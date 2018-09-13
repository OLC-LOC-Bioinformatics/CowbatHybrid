#!/usr/bin/env python

import os
import csv
import glob
from Bio import SeqIO


def create_combinedmetadata_report(assemblies_dir, reports_directory, metadata):
    if not os.path.isdir(reports_directory):
        os.makedirs(reports_directory)
    # Create basic assembly stats (N50, numcontigs, total length)
    basic_stats_report(assemblies_dir=assemblies_dir,
                       reports_directory=reports_directory)
    with open(os.path.join(reports_directory, 'combinedMetadata.csv'), 'w') as f:
        f.write('SampleName,N50,NumContigs,TotalLength,MeanInsertSize,AverageCoverageDepth,ReferenceGenome,RefGenomeAlleleMatches,16sPhylogeny,rMLSTsequenceType,MLSTsequencetype,MLSTmatches,coreGenome,SeroType,geneSeekrProfile,vtyperProfile,percentGC,TotalPredictedGenes,predictedgenesover3000bp,predictedgenesover1000bp,predictedgenesover500bp,predictedgenesunder500bp,SequencingDate,Investigator,TotalClustersinRun,NumberofClustersPF,PercentOfClusters,LengthofForwardRead,LengthofReverseRead,Project,PipelineVersion\n')
        for sample in metadata.runmetadata.samples:
            # Find n50, num_contigs and totallength by parsing the basic stats report.
            with open(os.path.join(reports_directory, 'basic_stats.csv')) as csvfile:
                sample_found = False
                reader = csv.DictReader(csvfile)
                for row in reader:
                    if sample.name == row['SampleName']:
                        sample_found = True
                        n50 = row['N50']
                        num_contigs = row['NumContigs']
                        total_length = row['TotalLength']
                if sample_found is False:
                    n50 = 'ND'
                    num_contigs = 'ND'
                    total_length = 'ND'
            f.write('{},{},{},{},'.format(sample.name, n50, num_contigs, total_length))
            for i in range(13):
                f.write(',')
            f.write('{},{},{},{},{},'.format(sample.prodigal.predictedgenestotal,
                                             sample.prodigal.predictedgenesover3000bp,
                                             sample.prodigal.predictedgenesover1000bp,
                                             sample.prodigal.predictedgenesover500bp,
                                             sample.prodigal.predictedgenesunder500bp))
            for i in range(9):
                f.write(',')
            f.write('\n')


def basic_stats_report(assemblies_dir, reports_directory):
    assembly_files = glob.glob(os.path.join(assemblies_dir, '*.fasta'))
    with open(os.path.join(reports_directory, 'basic_stats.csv'), 'w') as f:
        f.write('SampleName,N50,NumContigs,TotalLength\n')
    for assembly in assembly_files:
        total_length = 0
        contig_sizes = list()
        for contig in SeqIO.parse(assembly, 'fasta'):
            contig_length = len(contig)
            contig_sizes.append(contig_length)
            total_length += contig_length
        contig_sizes = sorted(contig_sizes, reverse=True)
        num_contigs = len(contig_sizes)
        length_so_far = 0
        n50 = 0
        i = 0
        while length_so_far <= (total_length * 0.5) and i < len(contig_sizes):
            length_so_far += contig_sizes[i]
            n50 = contig_sizes[i]
            i += 1
        with open(os.path.join(reports_directory, 'basic_stats.csv'), 'a+') as f:
            sample_name = os.path.split(assembly)[1].replace('.fasta', '')
            f.write('{},{},{},{}\n'.format(sample_name, n50, num_contigs, total_length))



