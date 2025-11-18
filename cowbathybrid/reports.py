#!/usr/bin/env python3

"""
Create a combined metadata report from the legacy combinedMetadata.csv file.
"""

# Standard imports
import csv
import os
from typing import Tuple

class Sample:
    def __init__(self, name, datastore, out_dir):
        self.name = name
        self.datastore = datastore
        self.general = General(out_dir)

class General:
    def __init__(self, out_dir):
        self.out_dir = out_dir

class RunMetadata:
    def __init__(self, samples):
        self.samples = samples

class Metadata:
    def __init__(self, runmetadata):
        self.runmetadata = runmetadata

def parse_flye_log(flye_log_path):
    """
    Parses the flye.log file to extract the Mean Coverage value.
    :param flye_log_path: Full path to the flye.log file
    :return: Mean Coverage value as a string
    """
    mean_coverage = 'ND'
    with open(flye_log_path, 'r') as log_file:
        for line in log_file:
            if 'Mean coverage:' in line:
                mean_coverage = line.split(':')[1].strip()
                break
    return mean_coverage

def parse_nanoplot_stats(nanoplot_stats_path) -> Tuple[str, str]:
    """
    Parses the NanoStats.txt file to extract the Mean read length value.
    :param nanoplot_stats_path: Full path to the NanoStats.txt file
    :return: Mean read length value as a string
    """
    mean_read_length = 'ND'
    mean_read_quality = 'ND'
    with open(nanoplot_stats_path, 'r') as stats_file:
        for line in stats_file:
            if 'Mean read length:' in line:
                mean_read_length = line.split(':')[1].strip()
            elif 'Mean read quality:' in line:
                mean_read_quality = line.split(':')[1].strip()
    return mean_read_length, mean_read_quality

def create_combinedmetadata_report(assemblies_dir, reports_directory, metadata):
    if not os.path.isdir(reports_directory):
        os.makedirs(reports_directory)

    legacy_metadata_path = os.path.join(reports_directory, 'legacy_combinedMetadata.csv')
    with open(legacy_metadata_path, 'r', encoding='utf-8') as combined_legacy:
        reader = csv.DictReader(combined_legacy)
        legacy_metadata = list(reader)

    for sample in metadata.runmetadata.samples:
        flye_log = os.path.join(sample.general.out_dir, 'flye', 'flye.log')
        nanoplot_stats = os.path.join(sample.general.out_dir, 'nanoplot', 'NanoStats.txt')

        mean_coverage = parse_flye_log(flye_log)
        mean_read_length, mean_read_quality = parse_nanoplot_stats(nanoplot_stats)

        for row in legacy_metadata:
            if row['SeqID'] == sample.name:
                row['Nanopore_coverage'] = mean_coverage
                row['Mean_read_length_nanopore'] = mean_read_length
                row['Mean_read_quality_nanopore'] = mean_read_quality

    with open(legacy_metadata_path, 'w', encoding='utf-8', newline='') as combined_legacy:
        fieldnames = legacy_metadata[0].keys()
        writer = csv.DictWriter(combined_legacy, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(legacy_metadata)