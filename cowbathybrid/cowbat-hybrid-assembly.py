#!/usr/bin/env python

from cowbathybrid.parsers import parse_hybrid_csv
from cowbathybrid.quality import run_nanoplot
from cowbathybrid import assemble
import multiprocessing
import argparse
import os


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
    args = parser.parse_args()

    # Parse the input CSV file we were given. This returns a list of SequenceFileInfo objects,
    # which have illumina_r1, illumina_r2, minion_reads, and outname as attributes.
    sequence_file_info_list = parse_hybrid_csv(args.input_csv)
    # Run NanoPlot on each of the MinION fastq files.
    # This creates a output_directory/samplename/nanoplot
    for sequence_file_info in sequence_file_info_list:
        nanoplot_outdir = os.path.join(args.output_directory, sequence_file_info.outname, 'nanoplot')
        if not os.path.isdir(nanoplot_outdir):
            os.makedirs(nanoplot_outdir)
        run_nanoplot(fastq_file=sequence_file_info.minion_reads,
                     threads=args.threads,
                     output_directory=nanoplot_outdir)

    # Give the list of sequence_file_info to the run_hybrid_assembly method.
    # It will run unicycler. Once this is done we should have a folder for each sample called outname in our
    # output directory. Within that, there will be a folder called unicycler that contains assembly.fasta
    assemble.run_hybrid_assembly(sequence_file_info_list=sequence_file_info_list,
                                 output_directory=args.output_directory,
                                 threads=args.threads)

    # Now we need to set up our metadata objects to do MLST/rMLST calling