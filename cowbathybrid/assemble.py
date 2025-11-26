#!/usr/bin/env python

import os
import gzip
import shutil
import logging
from cowbathybrid.command_runner import run_cmd


# Heng Li readfq - super fast!
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs);  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def trim_illumina(forward_reads, reverse_reads, output_directory, threads, logfile=None):
    forward_trimmed = os.path.join(output_directory, os.path.split(forward_reads.replace('.fastq.gz', '_trimmed.fastq.gz'))[1])
    reverse_trimmed = os.path.join(output_directory, os.path.split(reverse_reads.replace('.fastq.gz', '_trimmed.fastq.gz'))[1])
    cmd = 'bbduk.sh in={forward_reads} in2={reverse_reads} out={forward_trimmed} out2={reverse_trimmed} ' \
          'qtrim=w trimq=10 ref=adapters minlength=50 threads={threads}'.format(forward_reads=forward_reads,
                                                                                reverse_reads=reverse_reads,
                                                                                forward_trimmed=forward_trimmed,
                                                                                reverse_trimmed=reverse_trimmed,
                                                                                threads=threads)
    run_cmd(cmd, logfile=logfile)
    return forward_trimmed, reverse_trimmed


def correct_illumina(forward_reads, reverse_reads, output_directory, threads, logfile=None):
    forward_corrected = os.path.join(output_directory, os.path.split(forward_reads.replace('.fastq.gz', '_corrected.fastq.gz'))[1])
    reverse_corrected = os.path.join(output_directory, os.path.split(reverse_reads.replace('.fastq.gz', '_corrected.fastq.gz'))[1])
    cmd = 'tadpole.sh in={forward_reads} in2={reverse_reads} out={forward_corrected} out2={reverse_corrected} ' \
          'mode=correct threads={threads}'.format(forward_reads=forward_reads,
                                                  reverse_reads=reverse_reads,
                                                  forward_corrected=forward_corrected,
                                                  reverse_corrected=reverse_corrected,
                                                  threads=threads)
    run_cmd(cmd, logfile=logfile)
    return forward_corrected, reverse_corrected


#Madhu introduced the fly assembler to produce long contigs for producing hybrid assembly ##########################
def run_flye(long_reads, output_directory, threads, logfile=None):
    cmd = 'flye --nano-raw {long_reads} -t {threads} --out-dir {output_directory}'.format(long_reads=long_reads,
                                                                                               output_directory=output_directory,
											       threads=threads) #AC removed a runmode=runmode which was not defined
    run_cmd(cmd, logfile=logfile)
#    return flye_contigs AC removed on 230302 because it is not defined.
      

def run_unicycler(forward_reads, reverse_reads, long_reads, flye_contigs, output_directory, threads, logfile=None, conservative=False):
    if conservative:
        runmode = "conservative"
    else:
        runmode ="normal"
    cmd = 'unicycler -1 {forward_reads} -2 {reverse_reads} -l {long_reads} -o {output_directory} -t {threads} ' \
          '--no_correct --min_fasta_length 2000 --existing_long_read_assembly {flye_contigs} --keep 0 --mode {runmode}'.format(forward_reads=forward_reads, #AC changed keep to 1 on 230310
                                                                 reverse_reads=reverse_reads,
                                                                 flye_contigs=flye_contigs,
                                                                 long_reads=long_reads,
                                                                 output_directory=output_directory,
                                                                 threads=threads, runmode=runmode)
    run_cmd(cmd, logfile=logfile)


def run_porechop(minion_reads, output_directory, threads, logfile=None):
    chopped_reads = os.path.join(output_directory, 'minION_chopped.fastq.gz')
    cmd = 'porechop -i {minion_reads} -o {chopped_reads} -t {threads}'.format(minion_reads=minion_reads,
                                                                              chopped_reads=chopped_reads,
                                                                              threads=threads)
    run_cmd(cmd, logfile=logfile)
    return chopped_reads


def subsample_via_filtlong(minion_reads, illumina_forward, illumina_reverse, output_directory, target_bases=250000000, logfile=None):
    logging.info('Targeting {} bases for read subsampling...'.format(target_bases))
    filtered_reads = os.path.join(output_directory, 'length_filtered_reads.fastq')
    cmd = 'filtlong -t {target_bases} -1 {illumina_forward} -2 {illumina_reverse} {minion_reads} > {filtered_reads}'.format(target_bases=target_bases, illumina_forward=illumina_forward, illumina_reverse=illumina_reverse, minion_reads=minion_reads, filtered_reads=filtered_reads)
    run_cmd(cmd, logfile=logfile)
    return filtered_reads


# Andrew's manual read subsampling function -- not ideal for smaller plasmids
def subsample_minion_reads(minion_reads, output_directory, target_bases=250000000):
    # Ideally, would use SeqIO.index, but that doesn't work on gzipped files. Instead, use Heng Li's readfq
    # to iterate through the fastq file, and create a dictionary of read names/lengths
    read_name_length_dict = dict()
    logging.info('Targeting {} bases for read subsampling...'.format(target_bases))
    if minion_reads.endswith('.gz'):
        for name, seq, qual in readfq(gzip.open(minion_reads, 'rt')):
            read_name_length_dict[name] = len(seq)
    else:
        for name, seq, qual in readfq(open(minion_reads)):
            read_name_length_dict[name] = len(seq)
    # Now we have a nice dictionary with read names as keys and read lengths as values. Sort the dictionary by value
    # so that the longest reads are first.
    # This gives a list of tuples where index 0 of tuple is read name, index 1 is read length
    sorted_read_name_length = sorted(read_name_length_dict.items(), key=lambda kv: kv[1], reverse=True)
    # Now we iterate through the list of tuples until we hit our target number of bases.
    # Keep a set of read names we'll want to write.
    reads_to_write = set()
    total_bases_written = 0
    shortest_read_length = 'NA'
    for read in sorted_read_name_length:
        reads_to_write.add(read[0])
        total_bases_written += read[1]
        shortest_read_length = read[1]
        if total_bases_written > target_bases:
            break
    logging.info('Shortest read grabbed after subsampling: {} base pairs'.format(shortest_read_length))
    # With the set of read names we want to keep done as a set, do another parse through the fastq file and write
    # reads we want
    filtered_reads = os.path.join(output_directory, 'length_filtered_reads.fastq')
    with open(filtered_reads, 'w') as f:  # This might be real slow. Look into buffering the file writes.
        if minion_reads.endswith('.gz'):
            for name, seq, qual in readfq(gzip.open(minion_reads, 'rt')):
                if name in reads_to_write:
                    header = '@' + name
                    f.write(header + '\n' + seq + '\n+\n' + qual + '\n')
        else:
            for name, seq, qual in readfq(open(minion_reads)):
                if name in reads_to_write:
                    header = '@' + name
                    f.write(header + '\n' + seq + '\n+\n' + qual + '\n')
    return filtered_reads

def modify_assembly_headers(assembly_file, output_directory):
    """
    Modifies the headers in the assembly file from >1, >2, etc., to >contig_1, >contig_2, etc.
    :param assembly_file: The path to the .fasta file.
    :param output_directory: Directory where the .fasta file is located.
    """
    file_path = assembly_file
    temp_file_path = file_path + '.tmp'

    with open(file_path, 'r') as infile, open(temp_file_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                new_header = line.replace('>', '>contig_', 1)
                outfile.write(new_header)
            else:
                outfile.write(line)

    # Replace the original file with the modified one
    os.replace(temp_file_path, file_path)

def run_hybrid_assembly(long_reads,flye_contigs, forward_short_reads, reverse_short_reads, assembly_file,gfa_file, output_directory, filter_reads=None, conservative=False, threads=1):#AC added gfa_file 230314
    """
    Trims and corrects Illumina reads using BBDuk, removes addapters from MinION reads, and then runs unicycler.
    :param long_reads: Path to minION reads - uncorrected.
    :param forward_short_reads: Path to forward illumina reads.
    :param reverse_short_reads: Path to reverse illumina reads.
    :param assembly_file: The name/path of the final assembly file you want created
    :param output_directory: Directory where all the work will be done.
    :param filter_reads: If None, nothing happens. If a number, will subsample longest reads in the fastq file in order
    to reach that number of bases.
    :param threads: Number of threads to use for analysis. Defaults to 1.
    """
    if os.path.isfile(assembly_file):
        logging.info('Assembly for {} already exists...'.format(assembly_file))
        return  # Don't bother re-assembling something that's already assembled
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    logging.info('Beginning assembly of {}...'.format(assembly_file))
    logging.info('Trimming Illumina eads...')
    logfile = os.path.join(output_directory, 'hybrid_assembly_log.txt')
    forward_trimmed, reverse_trimmed = trim_illumina(forward_reads=forward_short_reads,
                                                     reverse_reads=reverse_short_reads,
                                                     output_directory=output_directory,
                                                     threads=threads,
                                                     logfile=logfile)
    logging.info('Correcting Illumina reads...')
    forward_corrected, reverse_corrected = correct_illumina(forward_reads=forward_trimmed,
                                                            reverse_reads=reverse_trimmed,
                                                            output_directory=output_directory,
                                                            threads=threads,
                                                            logfile=logfile)
    if filter_reads:
        logging.info('Subsampling MinION reads to choose the longest ones...')
        filtered_reads = subsample_via_filtlong(minion_reads=long_reads,
                                                illumina_forward=forward_short_reads,
                                                illumina_reverse=reverse_short_reads,
                                                output_directory=output_directory,
                                                target_bases=filter_reads)
        logging.info('Chopping adapters from minION reads...')
        chopped_reads = run_porechop(minion_reads=filtered_reads,
                                     output_directory=output_directory,
                                     threads=threads,
                                     logfile=logfile)
    else:
        filtered_reads = None
        logging.info('Chopping adapters from minION reads...')
        chopped_reads = run_porechop(minion_reads=long_reads,
                                     output_directory=output_directory,
                                     threads=threads,
                                     logfile=logfile)
        logging.info('running flye assembler from long reads...')
        run_flye(long_reads=chopped_reads,
                 output_directory=os.path.join(output_directory, 'flye'),
                 threads=threads,
                 logfile=logfile)
        logging.info('flye complete!')
        logging.info('Running Unicycler - this will take a while!')
        run_unicycler(
            forward_reads=forward_corrected,
            reverse_reads=reverse_corrected,
            flye_contigs=flye_contigs,
            long_reads=chopped_reads,
            output_directory=os.path.join(output_directory, 'unicycler'),
            threads=threads,
            logfile=logfile,
            conservative=conservative
        )
    logging.info('Unicycler complete!')
    shutil.copy(src=os.path.join(output_directory, 'unicycler', 'assembly.fasta'), dst=assembly_file)

    # Modify the assembly headers
    modify_assembly_headers(assembly_file, output_directory)

    shutil.copy(src=os.path.join(output_directory, 'unicycler', 'assembly.gfa'), dst=gfa_file)

    # Also remove the trimmed, corrected, and chopped files - they aren't necessary.
    os.remove(forward_trimmed)
    os.remove(forward_corrected)
    os.remove(reverse_corrected)
    os.remove(reverse_trimmed)
    os.remove(chopped_reads)
    if filter_reads:
        os.remove(filtered_reads)
