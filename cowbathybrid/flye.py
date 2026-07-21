#!/usr/bin/env python

from cowbathybrid.command_runner import run_cmd
import os
import logging


def run_flye(fastq_file, output_directory, threads, genome_size=None, asm_coverage=None):

    assembly_output = os.path.join(output_directory, 'assembly.fasta')

    if os.path.isfile(assembly_output):
        logging.info('Flye assembly already exists at %s, skipping...', assembly_output)
        return

    if (not os.path.isfile(fastq_file)) or os.path.getsize(fastq_file) == 0:
        raise FileNotFoundError(f"Flye input missing/empty: {fastq_file}")

    cmd = 'flye --nano-hq {fastq_file} -t {threads} --out-dir {output_directory}'.format(
        threads=threads,
        output_directory=output_directory,
        fastq_file=fastq_file
    )

    if genome_size:
        cmd += ' --genome-size {genome_size}'.format(genome_size=genome_size)

    if asm_coverage:
        cmd += ' --asm-coverage {asm_coverage}'.format(asm_coverage=asm_coverage)

    logging.info('Running Flye with input: %s', fastq_file)
    logging.info('Running Flye: %s', cmd)
    run_cmd(cmd)

    if not os.path.isfile(assembly_output):
        logging.error('Flye assembly failed - assembly.fasta not found at %s', assembly_output)
    else:
        logging.info('Flye assembly complete: %s', assembly_output)
