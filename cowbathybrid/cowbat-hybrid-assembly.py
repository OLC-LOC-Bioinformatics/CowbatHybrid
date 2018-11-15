#!/usr/bin/env python

from accessoryFunctions.accessoryFunctions import SetupLogging
from cowbathybrid.parsers import parse_hybrid_csv
from cowbathybrid.dependency_checks import check_dependencies
from cowbathybrid.quality import run_nanoplot
from cowbathybrid import assemble
from cowbat import assembly_typing
import multiprocessing
import argparse
import logging
import time
import os

# TODO: Adapter trimming - use PoreChop (now deprecated, so maybe not)? DownPore?

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly and perform some typing on hybrid MinION/Illumina data.')
    parser.add_argument('-i', '--input_csv',
                        required=True,
                        type=str,
                        help='Path to a CSV-formatted file with the following headers: MinION, Illumina_R1, '
                             'Illumina_R2, and OutName. For MinION, Illumina_R1, and Illumina_R2, the full '
                             'path to the read file for each sample should be present. OutName is what you want your '
                             'assembly to be called.')
    parser.add_argument('-r', '--referencefilepath',
                        required=True,
                        type=str,
                        help='Full path to folder containing reference databases.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help='Number of threads to run pipeline with. Defaults to number of cores on the system.')
    parser.add_argument('-o', '--output_directory',
                        type=str,
                        required=True,
                        help='Full path to directory where you want to store your outputs.')
    parser.add_argument('-verbose', '--verbose',
                        default=False,
                        action='store_true',
                        help='Activate this flag to get lots of debug output.')
    args = parser.parse_args()
    SetupLogging(debug=args.verbose)

    if check_dependencies() is False:
        quit(code=1)
    # Parse the input CSV file we were given. This returns a list of SequenceFileInfo objects,
    # which have illumina_r1, illumina_r2, minion_reads, and outname as attributes.
    sequence_file_info_list = parse_hybrid_csv(args.input_csv)
    # Run NanoPlot on each of the MinION fastq files.
    # This creates a output_directory/samplename/nanoplot
    for sequence_file_info in sequence_file_info_list:
        nanoplot_outdir = os.path.join(args.output_directory, sequence_file_info.outname, 'nanoplot')
        if not os.path.isdir(nanoplot_outdir):
            os.makedirs(nanoplot_outdir)
        logging.info('Running nanoplot on {}...'.format(sequence_file_info.outname))
        run_nanoplot(fastq_file=sequence_file_info.minion_reads,
                     threads=args.threads,
                     output_directory=nanoplot_outdir)

    best_assemblies_dir = os.path.join(args.output_directory, 'BestAssemblies')
    if not os.path.isdir(best_assemblies_dir):
        os.makedirs(best_assemblies_dir)
    for sequence_file_info in sequence_file_info_list:
        assemble.run_hybrid_assembly(long_reads=sequence_file_info.minion_reads,
                                     forward_short_reads=sequence_file_info.illumina_r1,
                                     reverse_short_reads=sequence_file_info.illumina_r2,
                                     output_directory=os.path.join(args.output_directory, sequence_file_info.outname, 'assembly'),
                                     threads=args.threads,
                                     assembly_file=os.path.join(best_assemblies_dir, sequence_file_info.outname + '.fasta'))

    # Much smarter way to do this - import Adam's assembly typing thingy and just use that.
    homepath = os.path.split(os.path.abspath(__file__))[0]   # No idea why this is necessary.
    typer = assembly_typing.Typing(start=time.time(),
                                   sequencepath=os.path.abspath(best_assemblies_dir),  # The CLARK portion of typing needs absolute path
                                   referencefilepath=args.referencefilepath,
                                   scriptpath=homepath)
    typer.main()
