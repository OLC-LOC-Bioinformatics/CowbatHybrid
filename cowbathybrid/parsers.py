#!/usr/bin/env python

import csv


class SequenceFileInfo:
    def __init__(self, minion_reads, illumina_r1, illumina_r2, outname):
        self.minion_reads = minion_reads
        self.illumina_r1 = illumina_r1
        self.illumina_r2 = illumina_r2
        self.outname = outname


def parse_hybrid_csv(csvfile):
    """
    Parses an input csv file so that the pipeline can know where each illumina file/minion file is and what to call
    each sample.
    :param csvfile: Full path to CSV file with headers MinION, Illumina_R1, Illumina_R2, and OutName
    :return: list of SequenceFileInfo objects
    """
    # TODO: Add checks that CSV file has proper headers
    sequence_file_info = list()
    with open(csvfile) as input_csv:
        reader = csv.DictReader(input_csv)
        for row in reader:
            sequence_file_info.append(SequenceFileInfo(minion_reads=row['MinION'],
                                                       illumina_r1=row['Illumina_R1'],
                                                       illumina_r2=row['Illumina_R2'],
                                                       outname=row['OutName']))
    return sequence_file_info
