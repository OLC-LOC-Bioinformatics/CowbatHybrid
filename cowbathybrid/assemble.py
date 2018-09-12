#!/usr/bin/env python

import os
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
