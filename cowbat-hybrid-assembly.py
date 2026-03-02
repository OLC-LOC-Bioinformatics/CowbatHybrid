#!/mnt/nas2/virtual_environments/cowbat_hybrid/bin/python3.9

"""
Wrapper script for the cowbat-hybrid pipeline.   This script will take a CSV file with the following headers:
MinION, Illumina_R1, Illumina_R2, and OutName.   For MinION, Illumina_R1, and Illumina_R2, the full path to the read
file for each sample should be present.  OutName is what you want your assembly to be called.  This script will run
nanoplot on the MinION reads, run flye on the MinION reads, and then run hybrid assemblies with Unicycler.   Finally,
it will run assembly typing on the assemblies with the cowbat pipeline.  
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
from cowbathybrid.command_runner import run_cmd
from cowbathybrid.dependency_checks import check_dependencies
from cowbathybrid.flye import run_flye
from cowbathybrid.parsers import parse_hybrid_csv
from cowbathybrid.quality import run_nanoplot
from cowbathybrid.version import __version__
from cowbathybrid import assemble
from cowbathybrid.reports import Metadata, create_combinedmetadata_report, Sample, RunMetadata

__author__ = 'Mathu Malar'

if __name__ == '__main__':
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Assembly and typing on hybrid MinION/Illumina data.')
    parser.add_argument(
        '-i', '--input_csv',
        required=True,
        type=str,
        help='Path to a CSV-formatted file with the following headers: MinION, Illumina_R1, Illumina_R2, and OutName. '
             'For MinION, Illumina_R1, and Illumina_R2, the full path to the read file for each sample should be '
             'present. OutName is what you want your assembly to be called.'
    )
    parser.add_argument(
        '-r', '--referencefilepath',
        required=True,
        type=str,
        help='Full path to folder containing reference databases.'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=multiprocessing.cpu_count(),
        help='Number of threads to run pipeline with. Defaults to number of cores on the system.'
    )
    parser.add_argument(
        '-o', '--output_directory',
        type=str,
        required=True,
        help='Full path to directory where you want to store your outputs.'
    )
    parser.add_argument(
        '-verbose', '--verbose',
        default=False,
        action='store_true',
        help='Activate this flag to get lots of debug output.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=__version__
    )
    parser.add_argument(
        '--asm-coverage',
        required=False,
        type=int,
        default=None,               # FIXED: None instead of False
        dest='asm_coverage',        # FIXED: explicit dest so args.asm_coverage works
        help='Reduced coverage for initial disjointig assembly (e.g. 50). '
             'Optional - set if flye is not producing any disjointigs. '
             'Recommended to use together with --genome-size.'
    )
    parser.add_argument(
        '-g', '--genome-size',
        required=False,
        type=str,                   # FIXED: str instead of int (values like 5m, 1.7m, 500k)
        default=None,               # FIXED: None instead of False
        dest='genome_size',         # FIXED: explicit dest so args.genome_size works
        help='For flye - Estimated genome size of assembled organism '
             '(e.g. 5m, 3m, 1.7m, 500k). '
             'Optional - if not provided, Flye estimates automatically. '
             'Required if using --asm-coverage flag.'
    )
    parser.add_argument(
        '-f', '--filter_reads',
        type=int,
        default=None,
        help='Unicycler can be pretty darn slow if given very large read sets. With this option, you can specify the '
             'number of long read bases you want to use. Aiming for 40X-50X depth for your target organism seems to '
             'work pretty well.'
    )
    parser.add_argument(
        '-c', '--conservative',
        default=False,
        action="store_true",
        help="Run Unicycler in conservative mode. Get more contigs, but fewer mis-assemblies."
    )
    args = parser.parse_args()
    SetupLogging(debug=args.verbose)

    if check_dependencies() is False:
        quit(code=1)

    # Log genome size and asm coverage if provided
    if args.genome_size:
        logging.info('Genome size set to: %s', args.genome_size)
    if args.asm_coverage:
        logging.info('Flye assembly coverage set to: %s', args.asm_coverage)

    # Parse the input CSV file we were given. This returns a list of SequenceFileInfo objects,
    # which have illumina_r1, illumina_r2, minion_reads, and outname as attributes.
    sequence_file_info_list = parse_hybrid_csv(args.input_csv)

    # Run NanoPlot on each of the MinION fastq files.
    # This creates an output_directory/samplename/nanoplot
    for sequence_file_info in sequence_file_info_list:
        nanoplot_out_dir = os.path.join(args.output_directory, sequence_file_info.outname, 'nanoplot')
        if not os.path.isdir(nanoplot_out_dir):
            os.makedirs(nanoplot_out_dir)
        logging.info('Running nanoplot on %s...', sequence_file_info.outname)
        run_nanoplot(
            fastq_file=sequence_file_info.minion_reads,
            threads=args.threads,
            output_directory=nanoplot_out_dir
        )

    # Run Flye assemblies to generate long contigs
    # Store flye output paths for each sample
    flye_outputs = {}
    for sequence_file_info in sequence_file_info_list:
        flye_out_dir = os.path.join(args.output_directory, sequence_file_info.outname, 'flye')
        if not os.path.isdir(flye_out_dir):
            os.makedirs(flye_out_dir)
        logging.info('Running flye on %s...', sequence_file_info.outname)
        run_flye(
            fastq_file=sequence_file_info.minion_reads,
            threads=args.threads,
            output_directory=flye_out_dir,
            genome_size=args.genome_size,       # FIXED: uncommented and corrected
            asm_coverage=args.asm_coverage      # FIXED: uncommented and corrected
        )
        # Store the flye output path for this specific sample
        flye_outputs[sequence_file_info.outname] = os.path.join(flye_out_dir, 'assembly.fasta')

    # Now run assemblies - intermediate files will be in output_directory/samplename/assembly
    # completed fasta file will be in output_directory/samplename.fasta
    best_assemblies_dir = os.path.join(args.output_directory, 'BestAssemblies')
    gfa_files_dir = os.path.join(args.output_directory, 'GFA_files')

    # Create the output directories if they don't exist
    os.makedirs(best_assemblies_dir, exist_ok=True)
    os.makedirs(gfa_files_dir, exist_ok=True)

    logging.info('Running Unicycler on samples...')

    # Run Unicycler on each of the samples with the correct flye output for that sample
    for sequence_file_info in sequence_file_info_list:
        assemble.run_hybrid_assembly(
            long_reads=sequence_file_info.minion_reads,
            flye_contigs=flye_outputs[sequence_file_info.outname],
            forward_short_reads=sequence_file_info.illumina_r1,
            reverse_short_reads=sequence_file_info.illumina_r2,
            output_directory=os.path.join(args.output_directory, sequence_file_info.outname, 'assembly'),
            threads=args.threads,
            assembly_file=os.path.join(best_assemblies_dir, sequence_file_info.outname + '.fasta'),
            gfa_file=os.path.join(gfa_files_dir, sequence_file_info.outname + '.gfa'),
            filter_reads=args.filter_reads,
            conservative=args.conservative
        )

    logging.info('Running assembly typing on samples...')
    home_path = os.path.split(os.path.abspath(__file__))[0]
    typer = assembly_typing.Typing(
        start=time.time(),
        sequencepath=os.path.abspath(best_assemblies_dir),
        referencefilepath=os.path.abspath(args.referencefilepath),
        scriptpath=home_path,
        debug=True
    )
    typer.main()

    # Clean up nested BestAssemblies directory created by COWBAT
    nested_best_assemblies = os.path.join(best_assemblies_dir, 'BestAssemblies')
    if os.path.isdir(nested_best_assemblies):
        logging.info('Cleaning up nested BestAssemblies directory...')
        shutil.rmtree(nested_best_assemblies)
        logging.info('Nested BestAssemblies removed.')

    logging.info('Creating reports...')
    samples = [
        Sample(
            name=sequence_file_info.outname,
            datastore='datastore',
            out_dir=os.path.join(args.output_directory, sequence_file_info.outname)
        )
        for sequence_file_info in sequence_file_info_list
    ]
    runmetadata = RunMetadata(samples=samples)
    metadata = Metadata(runmetadata=runmetadata)

    # Reports directory inside BestAssemblies/reports
    create_combinedmetadata_report(
        assemblies_dir=best_assemblies_dir,
        reports_directory=os.path.join(best_assemblies_dir, 'reports'),
        metadata=metadata
    )
    logging.info('Done!')