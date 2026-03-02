#!/usr/bin/env python

from cowbathybrid.command_runner import run_cmd
import os
import logging


# Ashley removed the --fastq_rich command, and changed flye to have --out-dir (was missing a -)
# Added optional genome_size and asm_coverage parameters (madhu 2026)
def run_flye(fastq_file, output_directory, threads, genome_size=None, asm_coverage=None):

    assembly_output = os.path.join(output_directory, 'assembly.fasta')

    # ── SKIP: if assembly.fasta already exists, flye has already run ──
    if os.path.isfile(assembly_output):
        logging.info('Flye assembly already exists at %s, skipping...', assembly_output)
        return
    # ──────────────────────────────────────────────────────────────────

    # Build base command
    cmd = 'flye --nano-hq {fastq_file} -t {threads} --out-dir {output_directory}'.format(
        threads=threads,
        output_directory=output_directory,
        fastq_file=fastq_file
    )

    # Only add --genome-size if provided (not all genomes are 5m)
    if genome_size:
        cmd += ' --genome-size {genome_size}'.format(genome_size=genome_size)

    # Only add --asm-coverage if provided
    if asm_coverage:
        cmd += ' --asm-coverage {asm_coverage}'.format(asm_coverage=asm_coverage)

    logging.info('Running Flye: %s', cmd)
    run_cmd(cmd)

    # Verify output was created
    if not os.path.isfile(assembly_output):
        logging.error('Flye assembly failed - assembly.fasta not found at %s', assembly_output)
    else:
        logging.info('Flye assembly complete: %s', assembly_output)