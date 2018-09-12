#!/usr/bin/env python

import os


def create_combinedmetadata_report(sequence_file_info_list, n50_dict, num_contigs_dict, total_length_dict,
                                   output_directory):
    reports_dir = os.path.join(output_directory, 'reports')
    if not os.path.isdir(reports_dir):
        os.makedirs(reports_dir)

    with open(os.path.join(reports_dir, 'combinedMetadata.csv'), 'w') as f:
        f.write('SampleName,N50,NumContigs,TotalLength,MeanInsertSize,AverageCoverageDepth,ReferenceGenome,RefGenomeAlleleMatches,16sPhylogeny,rMLSTsequenceType,MLSTsequencetype,MLSTmatches,coreGenome,SeroType,geneSeekrProfile,vtyperProfile,percentGC,TotalPredictedGenes,predictedgenesover3000bp,predictedgenesover1000bp,predictedgenesover500bp,predictedgenesunder500bp,SequencingDate,Investigator,TotalClustersinRun,NumberofClustersPF,PercentOfClusters,LengthofForwardRead,LengthofReverseRead,Project,PipelineVersion\n')
        for sequence_file_info in sequence_file_info_list:
            if sequence_file_info.outname in n50_dict:
                n50 = n50_dict[sequence_file_info.outname]
            else:
                n50 = 'ND'
            if sequence_file_info.outname in num_contigs_dict:
                num_contigs = num_contigs_dict[sequence_file_info.outname]
            else:
                num_contigs = 'ND'
            if sequence_file_info.outname in total_length_dict:
                total_length = total_length_dict[sequence_file_info.outname]
            else:
                total_length = 'ND'
            f.write('{},{},{},{},'.format(sequence_file_info.outname, n50, num_contigs, total_length))
            for i in range(27):
                f.write(',')
            f.write('\n')

