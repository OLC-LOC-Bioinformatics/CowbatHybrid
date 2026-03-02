#!/usr/bin/env python3

"""
Create a combined metadata report from the legacy combinedMetadata.csv file.
"""

# Standard imports
import csv
import os
import subprocess
import logging
from typing import Tuple, Dict

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
    Handles flye.log format:  "	Mean coverage:\t34"
    :param flye_log_path: Full path to the flye.log file
    :return: Mean Coverage value as a string
    """
    mean_coverage = 'ND'
    if not os.path.isfile(flye_log_path):
        logging.warning(f"Flye log not found: {flye_log_path}")
        return mean_coverage

    with open(flye_log_path, 'r') as log_file:
        for line in log_file:
            stripped = line.strip()
            # FIXED: startswith avoids matching timestamp colons
            # FIXED: split on 'Mean coverage:' not just ':'
            # Handles tab or space after colon
            if stripped.startswith('Mean coverage:'):
                try:
                    mean_coverage = stripped.split('Mean coverage:')[-1].strip()
                except (IndexError, ValueError):
                    mean_coverage = 'ND'
                break
    return mean_coverage

def parse_nanoplot_stats(nanoplot_stats_path) -> Tuple[str, str, str, str, str, str]:
    """
    Parses the NanoStats.txt file to extract nanopore read statistics.
    :param nanoplot_stats_path: Full path to the NanoStats.txt file
    :return: Tuple of (mean_read_length, mean_read_quality, median_read_length,
                       median_read_quality, number_of_reads, read_length_n50) as strings
    """
    mean_read_length    = 'ND'
    mean_read_quality   = 'ND'
    median_read_length  = 'ND'
    median_read_quality = 'ND'
    number_of_reads     = 'ND'
    read_length_n50     = 'ND'

    if not os.path.isfile(nanoplot_stats_path):
        logging.warning(f"NanoPlot stats not found: {nanoplot_stats_path}")
        return (mean_read_length, mean_read_quality,
                median_read_length, median_read_quality,
                number_of_reads, read_length_n50)

    with open(nanoplot_stats_path, 'r') as stats_file:
        for line in stats_file:
            stripped = line.strip()
            # Remove commas from numbers e.g. "3,585.6" -> "3585.6"
            if 'Mean read length:' in stripped:
                mean_read_length = stripped.split('Mean read length:')[-1].strip().replace(',', '')
            elif 'Mean read quality:' in stripped:
                mean_read_quality = stripped.split('Mean read quality:')[-1].strip().replace(',', '')
            elif 'Median read length:' in stripped:
                median_read_length = stripped.split('Median read length:')[-1].strip().replace(',', '')
            elif 'Median read quality:' in stripped:
                median_read_quality = stripped.split('Median read quality:')[-1].strip().replace(',', '')
            elif 'Number of reads:' in stripped:
                number_of_reads = stripped.split('Number of reads:')[-1].strip().replace(',', '')
            elif 'Read length N50:' in stripped:
                read_length_n50 = stripped.split('Read length N50:')[-1].strip().replace(',', '')

    return (mean_read_length, mean_read_quality,
            median_read_length, median_read_quality,
            number_of_reads, read_length_n50)

def run_quast(assembly_file, output_directory) -> Dict[str, str]:
    """
    Runs QUAST on the assembly file and extracts N50, NumContigs, TotalLength, and PercentGC.
    :param assembly_file: Full path to the assembly.fasta file
    :param output_directory: Directory where QUAST output will be stored
    :return: Dictionary with N50, NumContigs, TotalLength, and PercentGC values
    """
    quast_results = {
        'N50': 'ND',
        'NumContigs': 'ND',
        'TotalLength': 'ND',
        'PercentGC': 'ND'
    }

    if not os.path.isfile(assembly_file):
        logging.warning(f"Assembly file not found: {assembly_file}")
        return quast_results

    # Create QUAST output directory
    quast_out_dir = os.path.join(output_directory, 'quast')
    if not os.path.isdir(quast_out_dir):
        os.makedirs(quast_out_dir)

    # Run QUAST
    try:
        cmd = ['quast', assembly_file, '-o', quast_out_dir, '--silent']
        logging.info(f"Running QUAST: {' '.join(cmd)}")

        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Parse QUAST report
        report_file = os.path.join(quast_out_dir, 'report.tsv')
        if os.path.isfile(report_file):
            with open(report_file, 'r') as f:
                for line in f:
                    if line.startswith('# contigs'):
                        quast_results['NumContigs'] = line.split('\t')[1].strip()
                    elif line.startswith('Total length'):
                        quast_results['TotalLength'] = line.split('\t')[1].strip()
                    elif line.startswith('N50'):
                        quast_results['N50'] = line.split('\t')[1].strip()
                    elif line.startswith('GC (%)'):
                        quast_results['PercentGC'] = line.split('\t')[1].strip()
            logging.info(f"QUAST results for {assembly_file}: {quast_results}")
        else:
            logging.error(f"QUAST report not found: {report_file}")

    except subprocess.CalledProcessError as e:
        logging.error(f"QUAST failed for {assembly_file}: {e}")
        if hasattr(e, 'stdout') and e.stdout:
            logging.error(f"STDOUT: {e.stdout}")
        if hasattr(e, 'stderr') and e.stderr:
            logging.error(f"STDERR: {e.stderr}")
    except Exception as e:
        logging.error(f"Unexpected error running QUAST for {assembly_file}: {e}")

    return quast_results

def create_combinedmetadata_report(assemblies_dir, reports_directory, metadata):
    """
    Create combined metadata report with assembly metrics.

    :param assemblies_dir: Directory containing assembly FASTA files
    :param reports_directory: Directory where reports will be saved
    :param metadata: Metadata object containing sample information
    """
    if not os.path.isdir(reports_directory):
        os.makedirs(reports_directory)

    legacy_metadata_path = os.path.join(reports_directory, 'legacy_combinedMetadata.csv')

    # Canonical base fields - these must always be present in the output
    # New nanopore columns added at the end
    base_fields = [
        'SeqID',
        'Nanopore_coverage',
        'Mean_read_length_nanopore',
        'Mean_read_quality_nanopore',
        'N50',
        'NumContigs',
        'TotalLength',
        'PercentGC',
        'Median_read_length_nanopore',
        'Median_read_quality_nanopore',
        'Number_of_reads_nanopore',
        'Read_length_N50_nanopore',
    ]

    # Check if legacy file exists, create if not
    if os.path.isfile(legacy_metadata_path):
        with open(legacy_metadata_path, 'r', encoding='utf-8') as combined_legacy:
            reader = csv.DictReader(combined_legacy)
            legacy_metadata = list(reader)
            # Preserve the originally-parsed headers
            original_fieldnames = list(reader.fieldnames) if reader.fieldnames else list(base_fields)
        logging.info(f"Loaded existing metadata for {len(legacy_metadata)} samples")

        # Upgrade old legacy headers to include any new required fields
        # Handles legacy files created before new columns were added
        for f in base_fields:
            if f not in original_fieldnames:
                original_fieldnames.append(f)
                logging.info(f"Upgraded legacy header: added missing field '{f}'")

    else:
        logging.info("No legacy metadata file found, creating new file...")
        legacy_metadata = []
        original_fieldnames = list(base_fields)
        for sample in metadata.runmetadata.samples:
            legacy_metadata.append({
                'SeqID':                       sample.name,
                'Nanopore_coverage':           'ND',
                'Mean_read_length_nanopore':   'ND',
                'Mean_read_quality_nanopore':  'ND',
                'N50':                         'ND',
                'NumContigs':                  'ND',
                'TotalLength':                 'ND',
                'PercentGC':                   'ND',
                'Median_read_length_nanopore': 'ND',
                'Median_read_quality_nanopore':'ND',
                'Number_of_reads_nanopore':    'ND',
                'Read_length_N50_nanopore':    'ND',
            })

    # Create a set of existing sample names
    sample_names = {sample.name for sample in metadata.runmetadata.samples}

    # Create a list to hold updated metadata rows
    updated_metadata = []

    # Process each sample
    for sample in metadata.runmetadata.samples:
        flye_log      = os.path.join(sample.general.out_dir, 'flye', 'flye.log')
        nanoplot_stats = os.path.join(sample.general.out_dir, 'nanoplot', 'NanoStats.txt')
        assembly_file  = os.path.join(assemblies_dir, sample.name + '.fasta')

        logging.info(f"Processing metrics for sample: {sample.name}")

        mean_coverage = parse_flye_log(flye_log)

        (mean_read_length, mean_read_quality,
         median_read_length, median_read_quality,
         number_of_reads, read_length_n50) = parse_nanoplot_stats(nanoplot_stats)

        # Run QUAST on the assembly file
        quast_results = run_quast(assembly_file, sample.general.out_dir)

        # Update or add row
        found = False
        for row in legacy_metadata:
            # Only update rows for samples in our current run
            if row['SeqID'] not in sample_names:
                continue
            if row['SeqID'] == sample.name:
                row['Nanopore_coverage']            = mean_coverage
                row['Mean_read_length_nanopore']    = mean_read_length
                row['Mean_read_quality_nanopore']   = mean_read_quality
                row['N50']                          = quast_results['N50']
                row['NumContigs']                   = quast_results['NumContigs']
                row['TotalLength']                  = quast_results['TotalLength']
                row['PercentGC']                    = quast_results['PercentGC']
                row['Median_read_length_nanopore']  = median_read_length
                row['Median_read_quality_nanopore'] = median_read_quality
                row['Number_of_reads_nanopore']     = number_of_reads
                row['Read_length_N50_nanopore']     = read_length_n50
                found = True

                # Ensure the updated row has all fieldnames (fill missing with ND)
                for f in original_fieldnames:
                    row.setdefault(f, 'ND')

                updated_metadata.append(row)
                break

        if not found:
            logging.warning(f"Sample {sample.name} not found in existing metadata, adding new entry")
            # Build a new row that contains all the original headers (fill with 'ND')
            new_row = {f: 'ND' for f in original_fieldnames}
            new_row.update({
                'SeqID':                        sample.name,
                'Nanopore_coverage':            mean_coverage,
                'Mean_read_length_nanopore':    mean_read_length,
                'Mean_read_quality_nanopore':   mean_read_quality,
                'N50':                          quast_results['N50'],
                'NumContigs':                   quast_results['NumContigs'],
                'TotalLength':                  quast_results['TotalLength'],
                'PercentGC':                    quast_results['PercentGC'],
                'Median_read_length_nanopore':  median_read_length,
                'Median_read_quality_nanopore': median_read_quality,
                'Number_of_reads_nanopore':     number_of_reads,
                'Read_length_N50_nanopore':     read_length_n50,
            })
            updated_metadata.append(new_row)

    # Write updated metadata
    with open(legacy_metadata_path, 'w', encoding='utf-8', newline='') as combined_legacy:
        if updated_metadata:
            # Preserve original legacy file order: substitute updated rows where present,
            # then append any new samples (not in legacy) in the run's sample order
            updated_map = {u['SeqID']: u for u in updated_metadata}
            out_rows = []
            for row in legacy_metadata:
                seqid = row.get('SeqID')
                if seqid in updated_map:
                    out_rows.append(updated_map.pop(seqid))
                else:
                    out_rows.append(row)

            # Any remaining updated_map entries are new samples -> append in sample order
            for sample in metadata.runmetadata.samples:
                if sample.name in updated_map:
                    out_rows.append(updated_map.pop(sample.name))

            # Ensure every row has all fieldnames (fill missing with ND)
            for r in out_rows:
                for f in original_fieldnames:
                    r.setdefault(f, 'ND')

            # extrasaction='ignore' is a safety belt in case any row has unexpected keys
            writer = csv.DictWriter(combined_legacy, fieldnames=original_fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(out_rows)
            logging.info(f"Metadata report written to: {legacy_metadata_path}")
        else:
            logging.warning("No metadata to write")