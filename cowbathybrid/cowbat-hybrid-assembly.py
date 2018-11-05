#!/usr/bin/env python

from accessoryFunctions.accessoryFunctions import SetupLogging
from cowbathybrid.reports import create_combinedmetadata_report
from cowbathybrid.parsers import parse_hybrid_csv
from cowbathybrid.metadata_setup import Metadata
from cowbathybrid.dependency_checks import check_dependencies
from cowbathybrid.quality import run_nanoplot
from spadespipeline.typingclasses import Prophages, Univec
from spadespipeline.mobrecon import MobRecon
from spadespipeline.prodigal import Prodigal
from metagenomefilter import automateCLARK
import spadespipeline.sistr as sistr
from MASHsippr import mash as mash
import coreGenome.core as core
from cowbathybrid import assemble
from geneseekr.blast import BLAST
import multiprocessing
import argparse
import logging
import time
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly and perform some typing on hybrid MinION/Illumina data.')
    parser.add_argument('-i', '--input_csv',
                        required=True,
                        type=str,
                        help='Path to a CSV-formatted file with the following headers: MinION, Illumina_R1, '
                             'Illumina_R2, and OutName. For MinION, Illumina_R1, and Illumina_R2, the full '
                             'path to the read file for each sample should be present. OutName is what you want your '
                             'assembly to be called.')
    parser.add_argument('-r', '--referencefilepath',
                        required=True,
                        type=str,
                        help='Full path to folder containing reference databases.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help='Number of threads to run pipeline with. Defaults to number of cores on the system.')
    parser.add_argument('-o', '--output_directory',
                        type=str,
                        required=True,
                        help='Full path to directory where you want to store your outputs.')
    parser.add_argument('-verbose', '--verbose',
                        default=False,
                        action='store_true',
                        help='Activate this flag to get lots of debug output.')
    args = parser.parse_args()
    SetupLogging(debug=args.verbose)

    if check_dependencies() is False:
        quit(code=1)
    # Parse the input CSV file we were given. This returns a list of SequenceFileInfo objects,
    # which have illumina_r1, illumina_r2, minion_reads, and outname as attributes.
    sequence_file_info_list = parse_hybrid_csv(args.input_csv)
    # Run NanoPlot on each of the MinION fastq files.
    # This creates a output_directory/samplename/nanoplot
    for sequence_file_info in sequence_file_info_list:
        nanoplot_outdir = os.path.join(args.output_directory, sequence_file_info.outname, 'nanoplot')
        if not os.path.isdir(nanoplot_outdir):
            os.makedirs(nanoplot_outdir)
        logging.info('Running nanoplot on {}...'.format(sequence_file_info.outname))
        run_nanoplot(fastq_file=sequence_file_info.minion_reads,
                     threads=args.threads,
                     output_directory=nanoplot_outdir)

    best_assemblies_dir = os.path.join(args.output_directory, 'BestAssemblies')
    if not os.path.isdir(best_assemblies_dir):
        os.makedirs(best_assemblies_dir)
    for sequence_file_info in sequence_file_info_list:
        completed_assemblies = assemble.run_hybrid_assembly(long_reads=sequence_file_info.minion_reads,
                                                            forward_short_reads=sequence_file_info.illumina_r1,
                                                            reverse_short_reads=sequence_file_info.illumina_r2,
                                                            output_directory=os.path.join(args.output_directory, sequence_file_info.outname, 'assembly'),
                                                            threads=args.threads,
                                                            assembly_file=os.path.join(best_assemblies_dir, sequence_file_info.outname + '.fasta'))

    # All stuff for typing pipeline needs a 'MetadataObject', which keeps track of all the things and is what gets
    # passed to methods, modified, and then returned.
    metadata = Metadata(assemblies_dir=best_assemblies_dir,
                        starttime=time.time(),
                        logfile='log.txt',
                        outputdir=args.output_directory,
                        cpus=args.threads)
    metadata.strainer()
    metadata.reffilepath = args.referencefilepath
    metadata.reportpath = os.path.join(args.output_directory, 'reports')
    # Mash needs a trimmedcorrectedfastqfiles - hack that by making it equal to a list where only entry is the
    # best assembly file for that sample
    for sample in metadata.runmetadata.samples:
        sample.general.trimmedcorrectedfastqfiles = [sample.general.bestassemblyfile]
    # Try to determine what genus things are using MASH
    mash.Mash(inputobject=metadata,
              analysistype='mash')

    # Once we've mashed, we have a referencegenus attribute - other things expect a closestrefseqgenus attribute, so set
    # that too - they're pretty much the same thing.
    for sample in metadata.runmetadata.samples:
        sample.general.closestrefseqgenus = sample.general.referencegenus

    # Run all the assembly-based typing!
    Prodigal(metadata)

    # AMR
    metadata.targetpath = os.path.join(args.referencefilepath, 'resfinder')
    resfinder = BLAST(metadata, 'resfinder')
    resfinder.seekr()

    # Prophage
    metadata.targetpath = os.path.join(args.referencefilepath, 'prophages')
    prophages = Prophages(metadata, analysistype='prophages', cutoff=90)
    prophages.seekr()

    # Univec
    metadata.targetpath = os.path.join(args.referencefilepath, 'univec')
    univec = Univec(metadata, analysistype='univec', cutoff=80)
    univec.seekr()

    # Core genome
    metadata.targetpath = os.path.join(args.referencefilepath, 'coregenome')
    metadata.reffilepath = args.referencefilepath
    coregen = core.CoreGenome(metadata, analysistype='coregenome', genus_specific=True)
    coregen.seekr()
    core.AnnotatedCore(metadata)

    # Plasmids
    # MobRecon expects a resfinder_assembled attribute - this should be the same as resfinder, so set equal.
    for sample in metadata.runmetadata.samples:
        sample.resfinder_assembled = sample.resfinder

    mob = MobRecon(metadata=metadata.runmetadata.samples,
                   analysistype='mobrecon',
                   databasepath=metadata.reffilepath,
                   threads=args.threads,
                   logfile='log.txt',  # TODO: Set up actual logging.
                   reportpath=metadata.reportpath)
    mob.mob_recon()

    # Salmonella serotyping - genus should be taken care of by having run mash.
    sistr.Sistr(inputobject=metadata,
                analysistype='sistr')

    # CLARK
    automateCLARK.PipelineInit(inputobject=metadata,
                               extension='fasta',
                               light=True)

    # MLST

    # rMLST

    # Finally, create combined metadata report.
    create_combinedmetadata_report(assemblies_dir=os.path.join(args.output_directory, 'BestAssemblies'),
                                   reports_directory=os.path.join(args.output_directory, 'reports'),
                                   metadata=metadata)
    # reporter.Reporter(metadata)  TODO: Get the OLCTools version compatible with this?
