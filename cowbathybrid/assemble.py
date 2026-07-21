#!/usr/bin/env python

import os
import gzip
import shutil
import logging
from cowbathybrid.command_runner import run_cmd


def trim_illumina(forward_reads, reverse_reads, output_directory, threads, logfile=None):
    forward_trimmed = os.path.join(
        output_directory,
        os.path.split(forward_reads.replace('.fastq.gz', '_trimmed.fastq.gz'))[1]
    )
    reverse_trimmed = os.path.join(
        output_directory,
        os.path.split(reverse_reads.replace('.fastq.gz', '_trimmed.fastq.gz'))[1]
    )
    cmd = (
        f"bbduk.sh -Xmx120g -Xms60g in={forward_reads} in2={reverse_reads} "
        f"out={forward_trimmed} out2={reverse_trimmed} "
        f"qtrim=w trimq=10 ref=adapters minlength=50 threads={threads}"
    )
    run_cmd(cmd, logfile=logfile)
    return forward_trimmed, reverse_trimmed


def correct_illumina(forward_reads, reverse_reads, output_directory, threads, logfile=None):
    forward_corrected = os.path.join(
        output_directory,
        os.path.split(forward_reads.replace('.fastq.gz', '_corrected.fastq.gz'))[1]
    )
    reverse_corrected = os.path.join(
        output_directory,
        os.path.split(reverse_reads.replace('.fastq.gz', '_corrected.fastq.gz'))[1]
    )
    cmd = (
        f"tadpole.sh -Xmx180g in={forward_reads} in2={reverse_reads} "
        f"out={forward_corrected} out2={reverse_corrected} mode=correct threads={threads}"
    )
    run_cmd(cmd, logfile=logfile)
    return forward_corrected, reverse_corrected


def run_unicycler(forward_reads, reverse_reads, long_reads, flye_contigs, output_directory, threads,
                  logfile=None, conservative=False):
    runmode = "conservative" if conservative else "normal"
    cmd = (
        f"unicycler -1 {forward_reads} -2 {reverse_reads} -l {long_reads} "
        f"-o {output_directory} -t {threads} "
        f"--min_fasta_length 2000 "
        f"--existing_long_read_assembly {flye_contigs} --keep 0 --mode {runmode}"
    )
    run_cmd(cmd, logfile=logfile)


def run_porechop(minion_reads, output_directory, threads, logfile=None):
    chopped_reads = os.path.join(output_directory, 'minION_chopped.fastq.gz')
    if os.path.isfile(chopped_reads) and os.path.getsize(chopped_reads) > 0:
        logging.info("Using existing porechop output: %s", chopped_reads)
        return chopped_reads

    cmd = f"porechop -i {minion_reads} -o {chopped_reads} -t {threads}"
    run_cmd(cmd, logfile=logfile)

    if (not os.path.isfile(chopped_reads)) or os.path.getsize(chopped_reads) == 0:
        raise RuntimeError(f"Porechop output missing/empty: {chopped_reads}")

    return chopped_reads


def deduplicate_fastq_by_first_token(input_fastq_gz, output_fastq_gz):
    """
    Deduplicate FASTQ by first whitespace-delimited token in header (without leading '@').
    Keep first occurrence; drop subsequent duplicates.
    """
    seen = set()
    total = 0
    kept = 0
    dropped = 0

    with gzip.open(input_fastq_gz, 'rt') as fin, gzip.open(output_fastq_gz, 'wt') as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline()
            p = fin.readline()
            q = fin.readline()
            if not q:
                break

            total += 1
            key = h.strip().split()[0]
            if key.startswith('@'):
                key = key[1:]

            if key in seen:
                dropped += 1
                continue

            seen.add(key)
            kept += 1
            fout.write(h)
            fout.write(s)
            fout.write(p)
            fout.write(q)

    if (not os.path.isfile(output_fastq_gz)) or os.path.getsize(output_fastq_gz) == 0:
        raise RuntimeError(f"Dedup output missing/empty: {output_fastq_gz}")

    return total, kept, dropped


def qc_nanopore_reads_with_filtlong(minion_reads, output_directory, min_read_length=6000, keep_percent=95, logfile=None):
    dedup_reads = os.path.join(output_directory, 'minION_chopped_dedup.fastq.gz')
    total, kept, dropped = deduplicate_fastq_by_first_token(minion_reads, dedup_reads)
    logging.info("Dedup by first-token ID: total=%s, kept=%s, dropped_duplicates=%s", total, kept, dropped)

    tmp_fastq = os.path.join(output_directory, 'minION_chopped_qc.fastq')
    qc_reads = os.path.join(output_directory, 'minION_chopped_qc.fastq.gz')

    cmd = (
        f"filtlong --min_length {int(min_read_length)} "
        f"--keep_percent {float(keep_percent)} "
        f"{dedup_reads} > {tmp_fastq}"
    )
    run_cmd(cmd, logfile=logfile)

    if (not os.path.isfile(tmp_fastq)) or os.path.getsize(tmp_fastq) == 0:
        raise RuntimeError(
            f"filtlong produced empty output: {tmp_fastq}. "
            f"Try lower --min-read-length / higher --filtlong-keep-percent."
        )

    run_cmd(f"gzip -f {tmp_fastq}", logfile=logfile)

    if (not os.path.isfile(qc_reads)) or os.path.getsize(qc_reads) < 100:
        raise RuntimeError(f"Invalid filtlong gz output: {qc_reads}")

    return qc_reads, dedup_reads


def run_dnaapler_all(input_fasta, output_directory, threads=1, extra_args='', logfile=None):
    dnaapler_out_dir = os.path.join(output_directory, 'dnaapler')
    os.makedirs(dnaapler_out_dir, exist_ok=True)

    cmd = f"dnaapler all -i {input_fasta} -o {dnaapler_out_dir} -t {threads} {extra_args} --force"
    run_cmd(cmd, logfile=logfile)

    candidates = [
        os.path.join(dnaapler_out_dir, 'dnaapler_reoriented.fasta'),
        os.path.join(dnaapler_out_dir, 'reoriented.fasta'),
        os.path.join(dnaapler_out_dir, 'output.fasta')
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(f"Could not find dnaapler output fasta in {dnaapler_out_dir}")


def modify_assembly_headers(assembly_file):
    temp_file_path = assembly_file + '.tmp'
    with open(assembly_file, 'r') as infile, open(temp_file_path, 'w') as outfile:
        idx = 0
        for line in infile:
            if line.startswith('>'):
                idx += 1
                outfile.write(f'>contig_{idx}\n')
            else:
                outfile.write(line)
    os.replace(temp_file_path, assembly_file)


def run_hybrid_assembly(long_reads, flye_contigs, forward_short_reads, reverse_short_reads, assembly_file, gfa_file,
                        output_directory, filter_reads=None, conservative=False, threads=1,
                        min_read_length=6000, filtlong_keep_percent=95,
                        run_dnaapler=False, dnaapler_args=''):
    if os.path.isfile(assembly_file):
        logging.info('Assembly already exists: %s', assembly_file)
        return

    os.makedirs(output_directory, exist_ok=True)
    logfile = os.path.join(output_directory, 'hybrid_assembly_log.txt')

    forward_trimmed, reverse_trimmed = trim_illumina(
        forward_short_reads, reverse_short_reads, output_directory, threads, logfile
    )
    forward_corrected, reverse_corrected = correct_illumina(
        forward_trimmed, reverse_trimmed, output_directory, threads, logfile
    )

    # long_reads is already prepared in top-level workflow; don't re-run porechop/filtlong here
    long_reads_to_use = long_reads

    run_unicycler(
        forward_reads=forward_corrected,
        reverse_reads=reverse_corrected,
        long_reads=long_reads_to_use,
        flye_contigs=flye_contigs,
        output_directory=os.path.join(output_directory, 'unicycler'),
        threads=threads,
        logfile=logfile,
        conservative=conservative
    )

    unicycler_fasta = os.path.join(output_directory, 'unicycler', 'assembly.fasta')
    unicycler_gfa = os.path.join(output_directory, 'unicycler', 'assembly.gfa')

    if not os.path.isfile(unicycler_fasta):
        raise FileNotFoundError(
            f"Unicycler did not produce {unicycler_fasta}. "
            f"Check {logfile} and {os.path.join(output_directory, 'unicycler')}"
        )
    if not os.path.isfile(unicycler_gfa):
        raise FileNotFoundError(
            f"Unicycler did not produce {unicycler_gfa}. "
            f"Check {logfile} and {os.path.join(output_directory, 'unicycler')}"
        )

    shutil.copy(unicycler_fasta, assembly_file)
    shutil.copy(unicycler_gfa, gfa_file)

    if run_dnaapler:
        try:
            dnaapler_fasta = run_dnaapler_all(assembly_file, output_directory, threads, dnaapler_args, logfile)
            shutil.copy(dnaapler_fasta, assembly_file)
        except Exception as e:
            logging.warning("dnaapler failed (%s). Continuing with Unicycler assembly.", e)

    modify_assembly_headers(assembly_file)

    # cleanup
    for p in [forward_trimmed, reverse_trimmed, forward_corrected, reverse_corrected]:
        if os.path.isfile(p):
            os.remove(p)
