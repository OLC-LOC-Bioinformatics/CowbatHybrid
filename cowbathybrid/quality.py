#!/usr/bin/env python

from cowbathybrid.command_runner import run_cmd


def run_nanoplot(fastq_file, output_directory, threads):
    cmd = 'NanoPlot -t {threads} -o {output_directory} --fastq_rich {fastq_file}'.format(threads=threads,
                                                                                         output_directory=output_directory,
                                                                                         fastq_file=fastq_file)
    run_cmd(cmd)
