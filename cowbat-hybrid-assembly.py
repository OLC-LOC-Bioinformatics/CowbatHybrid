#!/usr/bin/env python

"""
Wrapper script for the cowbat-hybrid pipeline. This script will take a CSV file with the following headers:
MinION, Illumina_R1, Illumina_R2, and OutName. For MinION, Illumina_R1, and Illumina_R2, the full path to the read
file for each sample should be present. OutName is what you want your assembly to be called. This script will run
nanoplot on the MinION reads, run flye on the MinION reads, and then run hybrid assemblies with Unicycler. Finally,
it will run assembly typing on the assemblies with the cowbat pipeline.
"""

# Standard imports
import argparse
import logging
import multiprocessing
import os
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
        '-f', '--filter_reads',
        type=int,
        help='Unicycler can be pretty darn slow if given very large read sets. With this option, you can specify the '
             'number of long read bases you want to use. Aiming for 40X-50X depth for your target organism seems to '
             'work pretty well.'
    )
    parser.add_argument(
        '-c', '--conservative',
        default=False,
        action="store_true",
        help="Run Unicycler in conservative mode. Get more contigs, but fewer mis-assemblies"
    )
    args = parser.parse_args()
    SetupLogging(debug=args.verbose)

    if check_dependencies() is False:
        quit(code=1)
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

    # Initialize the output directory for the flye assemblies
    flye_out_dir = None

    ###Madhu: run fly assemblies to generate long contigs####
    for sequence_file_info in sequence_file_info_list:
        flye_out_dir = os.path.join(args.output_directory, sequence_file_info.outname, 'flye')
        if not os.path.isdir(flye_out_dir):
            os.makedirs(flye_out_dir)
        logging.info('Running flye on %s...', sequence_file_info.outname)
        run_flye(
            fastq_file=sequence_file_info.minion_reads,
            threads=args.threads,
            output_directory=flye_out_dir
        )

    # Now run assemblies - intermediate files will be in output_directory/samplename/assembly, completed fasta file
    # will be in output_directory/samplename.fasta
    best_assemblies_dir = os.path.join(args.output_directory)
    flye_output = os.path.join(flye_out_dir, 'assembly.fasta')
    gfa_files_dir = os.path.join(args.output_directory) #madhu added this

    # Create the output directory if it doesn't exist
    os.makedirs(best_assemblies_dir, exist_ok=True)

    logging.info('Running Unicycler on samples...')

    # Run Unicycler on each of the samples
    for sequence_file_info in sequence_file_info_list:
        assemble.run_hybrid_assembly(
            long_reads=sequence_file_info.minion_reads,
            flye_contigs=flye_output,
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

    logging.info('Creating reports...')
    samples = [
        Sample(name=sequence_file_info.outname, datastore='datastore', out_dir=os.path.join(args.output_directory, sequence_file_info.outname))
        for sequence_file_info in sequence_file_info_list
    ]
    runmetadata = RunMetadata(samples=samples)
    metadata = Metadata(runmetadata=runmetadata)
    create_combinedmetadata_report(
        assemblies_dir=best_assemblies_dir,
        reports_directory=os.path.join(os.path.abspath(best_assemblies_dir), 'reports'),
        metadata=metadata
    )
    logging.info('Done!')