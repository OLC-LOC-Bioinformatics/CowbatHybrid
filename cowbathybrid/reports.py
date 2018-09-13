#!/usr/bin/env python

import os
import glob
from Bio import SeqIO


def create_combinedmetadata_report(assemblies_dir, reports_directory, metadata):
    if not os.path.isdir(reports_directory):
        os.makedirs(reports_directory)
    # Create basic assembly stats (N50, numcontigs, total length)
    basic_stats_report(assemblies_dir=assemblies_dir,
                       reports_directory=reports_directory)
    # TODO: Actually write combined metadata report by parsing through all other reports.


def basic_stats_report(assemblies_dir, reports_directory):
    assembly_files = glob.glob(os.path.join(assemblies_dir, '*.fasta'))
    with open(os.path.join(reports_directory, 'basic_stats.csv'), 'w') as f:
        f.write('SampleName,N50,NumContigs,TotalLength\n')
    for assembly in assembly_files:
        total_length = 0
        contig_sizes = list()
        for contig in SeqIO.parse(assembly, 'fasta'):
            contig_length = len(contig)
            contig_sizes.append(contig_length)
            total_length += contig_length
        contig_sizes = sorted(contig_sizes, reverse=True)
        num_contigs = len(contig_sizes)
        length_so_far = 0
        n50 = 0
        i = 0
        while length_so_far <= (total_length * 0.5) and i < len(contig_sizes):
            length_so_far += contig_sizes[i]
            n50 = contig_sizes[i]
            i += 1
        with open(os.path.join(reports_directory, 'basic_stats.csv'), 'a+') as f:
            sample_name = os.path.split(assembly)[1].replace('.fasta', '')
            f.write('{},{},{},{}\n'.format(sample_name, n50, num_contigs, total_length))



