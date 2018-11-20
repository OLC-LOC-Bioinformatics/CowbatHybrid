#!/usr/bin/env python

import os
import shutil
import logging
from cowbathybrid.command_runner import run_cmd


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


def run_porechop(minion_reads, output_directory, threads):
    chopped_reads = os.path.join(output_directory, 'minION_chopped.fastq.gz')
    cmd = 'porechop -i {minion_reads} -o {chopped_reads} -t {threads}'.format(minion_reads=minion_reads,
                                                                              chopped_reads=chopped_reads,
                                                                              threads=threads)
    run_cmd(cmd)
    return chopped_reads


def run_hybrid_assembly(long_reads, forward_short_reads, reverse_short_reads, assembly_file, output_directory, threads):
    """
    Trims and corrects Illumina reads using BBDuk, removes addapters from MinION reads, and then runs unicycler.
    :param long_reads: Path to minION reads - uncorrected.
    :param forward_short_reads: Path to forward illumina reads.
    :param reverse_short_reads: Path to reverse illumina reads.
    :param assembly_file: The name/path of the final assembly file you want created
    :param output_directory: Directory where all the work will be done.
    :param threads: Number of threads to use for analysis
    """
    if os.path.isfile(assembly_file):
        logging.info('Assembly for {} already exists...'.format(assembly_file))
        return  # Don't bother re-assembling something that's already assembled
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    logging.info('Beginning assembly of {}...'.format(assembly_file))
    logging.info('Trimming Illumina eads...')
    forward_trimmed, reverse_trimmed = trim_illumina(forward_reads=forward_short_reads,
                                                     reverse_reads=reverse_short_reads,
                                                     output_directory=output_directory,
                                                     threads=threads)
    logging.info('Correcting Illumina reads...')
    forward_corrected, reverse_corrected = correct_illumina(forward_reads=forward_trimmed,
                                                            reverse_reads=reverse_trimmed,
                                                            output_directory=output_directory,
                                                            threads=threads)
    logging.info('Chopping adapters from minION reads...')
    chopped_reads = run_porechop(minion_reads=long_reads,
                                 output_directory=output_directory,
                                 threads=threads)
    logging.info('Running Unicycler - this will take a while!')
    run_unicycler(forward_reads=forward_corrected,
                  reverse_reads=reverse_corrected,
                  long_reads=chopped_reads,
                  output_directory=os.path.join(output_directory, 'unicycler'),
                  threads=threads)
    logging.info('Unicycler complete!')
    shutil.copy(src=os.path.join(output_directory, 'unicycler', 'assembly.fasta'), dst=assembly_file)

