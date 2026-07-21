#!/mnt/nas2/virtual_environments/cowbat_hybrid/bin/python3.9

"""
Wrapper script for the cowbat-hybrid pipeline. This script takes a CSV file with headers:
MinION, Illumina_R1, Illumina_R2, OutName.
"""

# Standard imports
import argparse
import logging
import multiprocessing
import os
import shutil
import time

# Third-party imports
from cowbat import assembly_typing
from olctools.accessoryFunctions.accessoryFunctions import SetupLogging

# Local imports
from cowbathybrid.dependency_checks import check_dependencies
from cowbathybrid.flye import run_flye
from cowbathybrid.parsers import parse_hybrid_csv
from cowbathybrid.quality import run_nanoplot
from cowbathybrid.version import __version__
from cowbathybrid import assemble
from cowbathybrid.reports import Metadata, create_combinedmetadata_report, Sample, RunMetadata

__author__ = 'Mathu Malar'


def relocate_nested_reports(best_assemblies_dir, root_reports_dir):
    nested_reports = os.path.join(best_assemblies_dir, 'reports')
    if not os.path.isdir(nested_reports):
        return

    os.makedirs(root_reports_dir, exist_ok=True)
    logging.info("Relocating reports: %s -> %s", nested_reports, root_reports_dir)

    for item in os.listdir(nested_reports):
        src = os.path.join(nested_reports, item)
        dst = os.path.join(root_reports_dir, item)
        if os.path.exists(dst):
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            else:
                os.remove(dst)
        shutil.move(src, dst)

    shutil.rmtree(nested_reports)
    logging.info("Nested reports relocated.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly and typing on hybrid MinION/Illumina data.')
    parser.add_argument('-i', '--input_csv', required=True, type=str,
                        help='CSV with headers: MinION, Illumina_R1, Illumina_R2, OutName')
    parser.add_argument('-r', '--referencefilepath', required=True, type=str,
                        help='Full path to folder containing reference databases.')
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads. Defaults to all cores.')
    parser.add_argument('-o', '--output_directory', type=str, required=True,
                        help='Output directory.')
    parser.add_argument('-verbose', '--verbose', default=False, action='store_true',
                        help='Verbose logging.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('--asm-coverage', required=False, type=int, default=None, dest='asm_coverage')
    parser.add_argument('-g', '--genome-size', required=False, type=str, default=None, dest='genome_size')
    parser.add_argument('-f', '--filter_reads', type=int, default=None)
    parser.add_argument('-c', '--conservative', default=False, action='store_true')
    parser.add_argument('--min-read-length', type=int, default=6000, dest='min_read_length')
    parser.add_argument('--filtlong-keep-percent', type=float, default=95, dest='filtlong_keep_percent')
    parser.add_argument('--run-dnaapler', action='store_true', default=False, dest='run_dnaapler')
    parser.add_argument('--dnaapler-args', type=str, default='', dest='dnaapler_args')

    args = parser.parse_args()
    SetupLogging(debug=args.verbose)

    if check_dependencies() is False:
        raise SystemExit(1)

    sequence_file_info_list = parse_hybrid_csv(args.input_csv)

    best_assemblies_dir = os.path.join(args.output_directory, 'BestAssemblies')
    gfa_files_dir = os.path.join(args.output_directory, 'GFA_files')
    root_reports_dir = os.path.join(args.output_directory, 'reports')
    os.makedirs(best_assemblies_dir, exist_ok=True)
    os.makedirs(gfa_files_dir, exist_ok=True)
    os.makedirs(root_reports_dir, exist_ok=True)

    # 1) NanoPlot on raw reads
    for s in sequence_file_info_list:
        nanoplot_out_dir = os.path.join(args.output_directory, s.outname, 'nanoplot')
        os.makedirs(nanoplot_out_dir, exist_ok=True)
        logging.info('Running NanoPlot on %s...', s.outname)
        run_nanoplot(
            fastq_file=s.minion_reads,
            threads=args.threads,
            output_directory=nanoplot_out_dir
        )

    # 2) Prepare long reads once per sample (porechop + dedup + filtlong), then use for Flye + Unicycler
    flye_outputs = {}
    prepared_long_reads = {}

    for s in sequence_file_info_list:
        sample_dir = os.path.join(args.output_directory, s.outname)
        prep_dir = os.path.join(sample_dir, 'longread_prep')
        flye_out_dir = os.path.join(sample_dir, 'flye')
        os.makedirs(prep_dir, exist_ok=True)
        os.makedirs(flye_out_dir, exist_ok=True)

        prep_log = os.path.join(prep_dir, 'longread_prep.log')
        logging.info('Preparing Nanopore reads for %s (porechop + dedup + filtlong)...', s.outname)

        chopped_reads = assemble.run_porechop(
            minion_reads=s.minion_reads,
            output_directory=prep_dir,
            threads=args.threads,
            logfile=prep_log
        )

        # Always create/check deduped reads (first-token header dedup)
        dedup_reads = os.path.join(prep_dir, 'minION_chopped_dedup.fastq.gz')
        total, kept, dropped = assemble.deduplicate_fastq_by_first_token(chopped_reads, dedup_reads)
        logging.info(
            "Dedup stats for %s: total=%s kept=%s dropped_duplicates=%s",
            s.outname, total, kept, dropped
        )

        # Try filtlong; if it fails, use dedup reads (NOT raw chopped reads)
        try:
            qc_reads, _dedup_reads = assemble.qc_nanopore_reads_with_filtlong(
                minion_reads=chopped_reads,  # function dedups internally and then runs filtlong
                output_directory=prep_dir,
                min_read_length=args.min_read_length,
                keep_percent=args.filtlong_keep_percent,
                logfile=prep_log
            )
            long_reads_for_downstream = qc_reads
            logging.info('Using filtlong QC reads for %s: %s', s.outname, long_reads_for_downstream)
        except Exception as e:
            logging.warning("Filtlong failed for %s (%s). Falling back to deduplicated reads.", s.outname, e)
            long_reads_for_downstream = dedup_reads

        prepared_long_reads[s.outname] = long_reads_for_downstream

        # Flye now uses prepared reads (NOT raw input)
        logging.info('Running Flye on %s using prepared reads...', s.outname)
        run_flye(
            fastq_file=long_reads_for_downstream,
            threads=args.threads,
            output_directory=flye_out_dir,
            genome_size=args.genome_size,
            asm_coverage=args.asm_coverage
        )
        flye_outputs[s.outname] = os.path.join(flye_out_dir, 'assembly.fasta')

    # 3) Hybrid assembly
    logging.info('Running Unicycler hybrid assembly...')
    for s in sequence_file_info_list:
        assemble.run_hybrid_assembly(
            long_reads=prepared_long_reads[s.outname],
            flye_contigs=flye_outputs[s.outname],
            forward_short_reads=s.illumina_r1,
            reverse_short_reads=s.illumina_r2,
            output_directory=os.path.join(args.output_directory, s.outname, 'assembly'),
            threads=args.threads,
            assembly_file=os.path.join(best_assemblies_dir, s.outname + '.fasta'),
            gfa_file=os.path.join(gfa_files_dir, s.outname + '.gfa'),
            filter_reads=args.filter_reads,
            conservative=args.conservative,
            min_read_length=args.min_read_length,
            filtlong_keep_percent=args.filtlong_keep_percent,
            run_dnaapler=args.run_dnaapler,
            dnaapler_args=args.dnaapler_args
        )

    # 4) Typing
    logging.info('Running assembly typing...')
    home_path = os.path.split(os.path.abspath(__file__))[0]
    typer = assembly_typing.Typing(
        start=time.time(),
        sequencepath=os.path.abspath(best_assemblies_dir),
        referencefilepath=os.path.abspath(args.referencefilepath),
        scriptpath=home_path,
        debug=True
    )
    typer.main()

    # Cleanup accidental nested BestAssemblies folder if created
    nested_best_assemblies = os.path.join(best_assemblies_dir, 'BestAssemblies')
    if os.path.isdir(nested_best_assemblies):
        shutil.rmtree(nested_best_assemblies)

    # Ensure reports end up at output_root/reports
    relocate_nested_reports(best_assemblies_dir, root_reports_dir)

    # 5) Combined metadata report
    logging.info('Creating combined metadata report...')
    samples = [
        Sample(name=s.outname, datastore='datastore', out_dir=os.path.join(args.output_directory, s.outname))
        for s in sequence_file_info_list
    ]
    metadata = Metadata(runmetadata=RunMetadata(samples=samples))

    create_combinedmetadata_report(
        assemblies_dir=best_assemblies_dir,
        reports_directory=root_reports_dir,
        metadata=metadata
    )

    logging.info('Done! Reports available at %s', root_reports_dir)
