#!/usr/bin/env python

import os
import shutil
import logging
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


def generate_and_sort_bam(forward_reads, reverse_reads, assembly, output_bam):
    cmd = 'bbmap.sh in={forward_reads} in2={reverse_reads} out={output_bam} ref={assembly} ' \
          'nodisk'.format(forward_reads=forward_reads,
                          reverse_reads=reverse_reads,
                          assembly=assembly,
                          output_bam=output_bam.replace('.bam', '_unsorted.bam'))
    run_cmd(cmd, logfile='log.txt')
    cmd = 'samtools sort -o {output_bam} {unsorted_bam}'.format(output_bam=output_bam,
                                                                unsorted_bam=output_bam.replace('.bam', '_unsorted.bam'))
    run_cmd(cmd, logfile='log.txt')
    cmd = 'samtools index {}'.format(output_bam)
    run_cmd(cmd, logfile='log.txt')


def run_pilon(draft_assembly, forward_reads, reverse_reads, output_assembly, output_dir=os.getcwd(), max_pilon_rounds=10):
    pilon_rounds = 1
    num_changes_made = 8888858
    assembly_to_use = draft_assembly
    while pilon_rounds < max_pilon_rounds and num_changes_made > 0:
        logging.debug('Begin pilon round {}'.format(pilon_rounds))
        pilon_dir = os.path.join(output_dir, 'pilon_{}'.format(pilon_rounds))
        if not os.path.isdir(pilon_dir):
            os.makedirs(pilon_dir)
        generate_and_sort_bam(forward_reads=forward_reads,
                              reverse_reads=reverse_reads,
                              assembly=assembly_to_use,
                              output_bam=os.path.join(pilon_dir, 'pilon.bam'))
        cmd = 'pilon --genome {assembly} --bam {bam} --outdir {pilon_dir} --changes'.format(assembly=assembly_to_use,
                                                                                            bam=os.path.join(pilon_dir, 'pilon.bam'),
                                                                                            pilon_dir=pilon_dir)
        run_cmd(cmd, logfile='log.txt')
        with open(os.path.join(pilon_dir, 'pilon.changes')) as f:
            num_changes_made = len(f.readlines())
        logging.debug('Pilon round {} complete. Made {} changes.'.format(pilon_rounds, num_changes_made))
        assembly_to_use = os.path.join(pilon_dir, 'pilon.fasta')
        pilon_rounds += 1
    shutil.copy(src=assembly_to_use, dst=os.path.join(output_dir, output_assembly))


def run_hybrid_assembly(sequence_file_info_list, output_directory, threads):
    # TODO: Rename contigs - the headers end up with _pilon added however many times pilon ran.
    for sequence_file_info in sequence_file_info_list:
        if os.path.isfile(os.path.join(output_directory, 'BestAssemblies', sequence_file_info.outname + '.fasta')):
            continue
        if not os.path.isdir(os.path.join(output_directory, sequence_file_info.outname, 'assembly')):
            os.makedirs(os.path.join(output_directory, sequence_file_info.outname, 'assembly'))
        run_wtdbg2(long_reads=sequence_file_info.minion_reads,
                   output_name=os.path.join(output_directory, sequence_file_info.outname, 'assembly', 'wtdbg2_assembly'),
                   threads=threads)
        run_pilon(draft_assembly=os.path.join(output_directory, sequence_file_info.outname, 'assembly', 'wtdbg2_assembly.fasta'),
                  forward_reads=sequence_file_info.illumina_r1,
                  reverse_reads=sequence_file_info.illumina_r2,
                  output_dir=os.path.join(output_directory, sequence_file_info.outname, 'assembly'),
                  output_assembly='assembly.fasta')
        # TODO: Add cleanup of all pilon folders - lots of big BAM files and whatnot that don't need to stick around
    completed_assemblies = move_assemblies(sequence_file_info_list, output_directory)
    return completed_assemblies


def move_assemblies(sequence_file_info_list, output_directory):
    completed_assemblies = list()
    best_assemblies_dir = os.path.join(output_directory, 'BestAssemblies')
    if not os.path.isdir(best_assemblies_dir):
        os.makedirs(best_assemblies_dir)
    for sequence_file_info in sequence_file_info_list:
        hybrid_assembly = os.path.join(output_directory, sequence_file_info.outname, 'assembly', 'assembly.fasta')
        if os.path.isfile(hybrid_assembly):
            shutil.copy(src=hybrid_assembly, dst=os.path.join(best_assemblies_dir, sequence_file_info.outname + '.fasta'))
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


