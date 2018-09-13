#!/usr/bin/env python

from cowbathybrid.command_runner import run_cmd
import os


def run_nanoplot(fastq_file, output_directory, threads):
    if os.path.isfile(os.path.join(output_directory, 'NanoPlot-report.html')):
        return
    cmd = 'NanoPlot -t {threads} -o {output_directory} --fastq_rich {fastq_file}'.format(threads=threads,
                                                                                         output_directory=output_directory,
                                                                                         fastq_file=fastq_file)
    run_cmd(cmd)
