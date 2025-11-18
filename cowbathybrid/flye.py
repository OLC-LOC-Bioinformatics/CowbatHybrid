#!/usr/bin/env python

from cowbathybrid.command_runner import run_cmd
import os


#def run_flye(fastq_file, output_directory, threads):
#    if os.path.isfile(os.path.join(output_directory, 'flye-report.html')):
#        return
#    cmd = 'flye --nano-raw {long_reads} --plasmids -t {threads} -out-dir {output_directory} --fastq_rich {fastq_file}'.format(threads=threads,
#                                                                                         output_directory=output_directory, 
#                                                                                         long_reads=long_reads, #AC added the long_reads def
#                                                                                         fastq_file=fastq_file)
#    run_cmd(cmd)


#Ashley removed the --fastq_rich command, and changed flye to have --out-dir (was missing a -)... it also did not like using long_reads for the call??? I dunno.
def run_flye(fastq_file, output_directory, threads):
    if os.path.isfile(os.path.join(output_directory, 'flye-report.html')):
        return
    cmd = 'flye --nano-raw {fastq_file} -t {threads} --out-dir {output_directory}'.format(threads=threads,
                                                                                         output_directory=output_directory, 
                                                                                         fastq_file=fastq_file)
    run_cmd(cmd)

