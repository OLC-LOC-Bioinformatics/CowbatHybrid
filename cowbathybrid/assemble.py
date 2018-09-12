#!/usr/bin/env python

import os
import shutil
from Bio import SeqIO
from cowbathybrid.command_runner import run_cmd


def run_hybrid_assembly(sequence_file_info_list, output_directory, threads):
    for sequence_file_info in sequence_file_info_list:
        if not os.path.isdir(os.path.join(output_directory, sequence_file_info.outname)):
            os.makedirs(os.path.join(output_directory, sequence_file_info.outname))
        forward_trimmed, reverse_trimmed = trim_illumina(forward_reads=sequence_file_info.illumina_r1,
                                                         reverse_reads=sequence_file_info.illumina_r2,
                                                         output_directory=os.path.join(output_directory, sequence_file_info.outname),
                                                         threads=threads)
        forward_corrected, reverse_corrected = correct_illumina(forward_reads=forward_trimmed,
                                                                reverse_reads=reverse_trimmed,
                                                                output_directory=os.path.join(output_directory, sequence_file_info.outname),
                                                                threads=threads)
        run_unicycler(forward_reads=forward_corrected,
                      reverse_reads=reverse_corrected,
                      long_reads=sequence_file_info.minion_reads,
                      output_directory=os.path.join(output_directory, sequence_file_info.outname, 'unicycler'),
                      threads=threads)
    completed_assemblies = move_assemblies(sequence_file_info_list, output_directory)
    return completed_assemblies


def move_assemblies(sequence_file_info_list, output_directory):
    completed_assemblies = list()
    best_assemblies_dir = os.path.join(output_directory, 'BestAssemblies')
    if not os.path.isdir(best_assemblies_dir):
        os.makedirs(best_assemblies_dir)
    for sequence_file_info in sequence_file_info_list:
        unicycler_assembly = os.path.join(output_directory, sequence_file_info.outname, 'unicycler', 'assembly.fasta')
        if os.path.isfile(unicycler_assembly):
            shutil.copy(src=unicycler_assembly, dst=os.path.join(best_assemblies_dir, sequence_file_info.outname + '.fasta'))
            completed_assemblies.append(os.path.join(best_assemblies_dir, sequence_file_info.outname + '.fasta'))
    return completed_assemblies


def trim_illumina(forward_reads, reverse_reads, output_directory, threads):
    forward_trimmed = os.path.join(output_directory, os.path.split(forward_reads.replace('.fastq.gz', '_trimmed.fastq.gz'))[1])
    reverse_trimmed = os.path.join(output_directory, os.path.split(reverse_reads.replace('.fastq.gz', '_trimmed.fastq.gz'))[1])
    cmd = 'bbduk.sh in={forward_reads} in2={reverse_reads} out={forward_trimmed} out2={reverse_trimmed} ' \
          'qtrim=w trimq=10 ref=adapters minlength=50 threads={threads}'.format(forward_reads=forward_reads,
                                                                                reverse_reads=reverse_reads,
                                                                                forward_trimmed=forward_trimmed,
                                                                                reverse_trimmed=reverse_trimmed,
                                                                                threads=threads)
    run_cmd(cmd)
    return forward_trimmed, reverse_trimmed


def correct_illumina(forward_reads, reverse_reads, output_directory, threads):
    forward_corrected = os.path.join(output_directory, os.path.split(forward_reads.replace('.fastq.gz', '_corrected.fastq.gz'))[1])
    reverse_corrected = os.path.join(output_directory, os.path.split(reverse_reads.replace('.fastq.gz', '_corrected.fastq.gz'))[1])
    cmd = 'tadpole.sh in={forward_reads} in2={reverse_reads} out={forward_corrected} out2={reverse_corrected} ' \
          'mode=correct threads={threads}'.format(forward_reads=forward_reads,
                                                  reverse_reads=reverse_reads,
                                                  forward_corrected=forward_corrected,
                                                  reverse_corrected=reverse_corrected,
                                                  threads=threads)
    run_cmd(cmd)
    return forward_corrected, reverse_corrected


def run_unicycler(forward_reads, reverse_reads, long_reads, output_directory, threads):
    cmd = 'unicycler -1 {forward_reads} -2 {reverse_reads} -l {long_reads} -o {output_directory} -t {threads} ' \
          '--no_correct --min_fasta_length 1000 --keep 0'.format(forward_reads=forward_reads,
                                                                 reverse_reads=reverse_reads,
                                                                 long_reads=long_reads,
                                                                 output_directory=output_directory,
                                                                 threads=threads)
    run_cmd(cmd)


def find_n50(assembly_files):
    """
    Returns N50 size.
    :param assembly_files: List of paths to assemblies.
    :return: Dictionary with N50 as value and path to assembly as key.
    """
    n50_dict = dict()
    for assembly in assembly_files:
        total_length = 0
        contig_sizes = list()
        for contig in SeqIO.parse(assembly, 'fasta'):
            contig_length = len(contig)
            contig_sizes.append(contig_length)
            total_length += contig_length
        contig_sizes = sorted(contig_sizes, reverse=True)
        length_so_far = 0
        n50 = 0
        i = 0
        while length_so_far <= (total_length * 0.5) and i < len(contig_sizes):
            length_so_far += contig_sizes[i]
            n50 = contig_sizes[i]
            i += 1
        n50_dict[os.path.split(assembly)[1].replace('.fasta', '')] = n50
    return n50_dict


def find_total_length(assembly_files):
    total_length_dict = dict()
    for assembly in assembly_files:
        total_length = 0
        for contig in SeqIO.parse(assembly, 'fasta'):
            total_length += len(contig)
        total_length_dict[os.path.split(assembly)[1].replace('.fasta', '')] = total_length
    return total_length_dict


def find_num_contigs(assembly_files):
    contigs_dict = dict()
    for assembly in assembly_files:
        num_contigs = 0
        for contig in SeqIO.parse(assembly, 'fasta'):
            num_contigs += 1
        contigs_dict[os.path.split(assembly)[1].replace('fasta', '')] = num_contigs
    return contigs_dict

