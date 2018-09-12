#!/usr/bin/env python

import subprocess


def run_hybrid_assembly(sequence_file_info_list):
    for sequence_file_info in sequence_file_info_list:
        print('Trim Illumina!')
        print('Correct Illumina!')
        print('Run Unicycler')
