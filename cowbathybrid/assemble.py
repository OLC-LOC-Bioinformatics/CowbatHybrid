#!/usr/bin/env python

import os
import glob
import logging
from Bio import SeqIO
from cowbathybrid.command_runner import run_cmd


def run_wtdbg2(long_reads, output_name, threads=4):
    cmd = 'wtdbg2 -t {threads} -i {long_reads} -fo {output_name}'.format(long_reads=long_reads, output_name=output_name,
                                                                         threads=threads)
    run_cmd(cmd, logfile='log.txt')
    # This works in wtdbg2 v2.1, but in v2.2 this changes to .ctg.lay.gz - will need to be aware of this when
    # the conda package gets updated.
    cmd = 'wtpoa-cns -t {threads} -i {output_name}.ctg.lay -fo {output_name}.fasta'.format(output_name=output_name,
                                                                                           threads=threads)
    run_cmd(cmd, logfile='log.txt')


def generate_and_sort_bam(forward_reads, reverse_reads, assembly, output_bam, threads=1):
    cmd = 'bbmap.sh in={forward_reads} in2={reverse_reads} out={output_bam} ref={assembly} ' \
          'nodisk threads={threads}'.format(forward_reads=forward_reads,
                                            reverse_reads=reverse_reads,
                                            assembly=assembly,
                                            output_bam=output_bam.replace('.bam', '_unsorted.bam'),
                                            threads=threads)
    run_cmd(cmd, logfile='log.txt')
    cmd = 'samtools sort -o {output_bam} {unsorted_bam}'.format(output_bam=output_bam,
                                                                unsorted_bam=output_bam.replace('.bam', '_unsorted.bam'))
    run_cmd(cmd, logfile='log.txt')
    cmd = 'samtools index {}'.format(output_bam)
    run_cmd(cmd, logfile='log.txt')


def run_pilon(draft_assembly, forward_reads, reverse_reads, output_assembly, threads=1, output_dir=os.getcwd(), max_pilon_rounds=10):
    pilon_rounds = 1
    num_changes_made = 8888858
    assembly_to_use = draft_assembly
    while pilon_rounds <= max_pilon_rounds and num_changes_made > 0:
        logging.debug('Begin pilon round {}'.format(pilon_rounds))
        pilon_dir = os.path.join(output_dir, 'pilon_{}'.format(pilon_rounds))
        if not os.path.isdir(pilon_dir):
            os.makedirs(pilon_dir)
        generate_and_sort_bam(forward_reads=forward_reads,
                              reverse_reads=reverse_reads,
                              assembly=assembly_to_use,
                              output_bam=os.path.join(pilon_dir, 'pilon.bam'),
                              threads=threads)
        cmd = 'pilon --genome {assembly} --bam {bam} --outdir {pilon_dir} --changes'.format(assembly=assembly_to_use,
                                                                                            bam=os.path.join(pilon_dir, 'pilon.bam'),
                                                                                            pilon_dir=pilon_dir)
        run_cmd(cmd, logfile='log.txt')
        with open(os.path.join(pilon_dir, 'pilon.changes')) as f:
            num_changes_made = len(f.readlines())
        logging.debug('Pilon round {} complete. Made {} changes.'.format(pilon_rounds, num_changes_made))
        assembly_to_use = os.path.join(pilon_dir, 'pilon.fasta')
        pilon_rounds += 1
    rename_contigs_and_copy_sequences(input_fasta=assembly_to_use, output_fasta=output_assembly)


def rename_contigs_and_copy_sequences(input_fasta, output_fasta):
    """
    Pilon adds a _pilon to each contig for each round done. That looks stupid, so rewrite the files.
    :param input_fasta: Path to FASTA that's been corrected lots of times by pilon.
    :param output_fasta: Path to output fasta file.
    """
    sequences_to_write = list()
    for sequence in SeqIO.parse(input_fasta, 'fasta'):
        sequence.id = sequence.id.split('_')[0]
        sequences_to_write.append(sequence)
    SeqIO.write(sequences=sequences_to_write, handle=output_fasta, format='fasta')


def run_hybrid_assembly(long_reads, forward_short_reads, reverse_short_reads, assembly_file, output_directory, threads,
                        keep_bams=False):
    """
    Runs an assembly using nanopore reads using wtdbg2, which claims to be pretty much as good as Canu, but runs in
    about a minute. Then polishes the assembly lots using pilon - either 10 rounds or zero changes, whichever comes
    first. May need to look at these and see if they're actually good settings.
    :param long_reads: Path to minION reads - uncorrected.
    :param forward_short_reads: Path to forward illumina reads.
    :param reverse_short_reads: Path to reverse illumina reads.
    :param assembly_file: The name/path of the final assembly file you want created
    :param output_directory: Directory where all the work will be done.
    :param threads: Number of threads to use for analysis
    :param keep_bams: Pilon needs a bamfile at each step. Set to False to have these bamfiles get deleted, or true
    to have the bamfiles kept in case you want to take a closer look at them.
    """
    # TODO: Rename contigs - the headers end up with _pilon added however many times pilon ran.
    if os.path.isfile(assembly_file):
        return  # Don't bother re-assembling something that's already assembled
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    run_wtdbg2(long_reads=long_reads,
               output_name=os.path.join(output_directory, 'wtdbg2_assembly'),
               threads=threads)
    run_pilon(draft_assembly=os.path.join(output_directory, 'wtdbg2_assembly.fasta'),
              forward_reads=forward_short_reads,
              reverse_reads=reverse_short_reads,
              output_dir=output_directory,
              output_assembly=assembly_file,
              threads=threads)
    # Clean up all the bam files you generate when running pilon, unless it's specified that they should be kept
    if keep_bams is False:
        bamfiles = glob.glob(os.path.join(output_directory, 'pilon_*', '*.bam*'))
        for bamfile in bamfiles:
            os.remove(bamfile)


