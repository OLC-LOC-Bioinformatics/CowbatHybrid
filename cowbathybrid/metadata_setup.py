#!/usr/bin/env python

from accessoryFunctions.accessoryFunctions import GenObject, MetadataObject
import glob
import os


class Metadata:

    def __init__(self, assemblies_dir, logfile, starttime, outputdir, cpus):
        self.assemblies_dir = assemblies_dir
        self.logfile = logfile
        self.starttime = starttime
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = list()
        self.outputdir = outputdir
        self.cpus = cpus
        self.targetpath = ''
        self.reportpath = ''

    def strainer(self):
        """
        Locate all the FASTA files in the supplied sequence path. Create basic metadata objects for
        each sample
        """
        assert os.path.isdir(self.assemblies_dir), 'Cannot locate sequence path as specified: {}' \
            .format(self.assemblies_dir)
        # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
        # begins with .fa
        strains = sorted(glob.glob(os.path.join(self.assemblies_dir, '*.fa*'.format(self.assemblies_dir))))
        # Populate the metadata object. This object will be populated to mirror the objects created in the
        # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
        # of a pipeline
        for sample in strains:
            # Create the object
            metadata = MetadataObject()
            # Set the base file name of the sequence. Just remove the file extension
            filename = os.path.splitext(os.path.split(sample)[1])[0]
            # Set the .name attribute to be the file name
            metadata.name = filename
            # Create the .general attribute
            metadata.general = GenObject()
            metadata.commands = GenObject()
            # Set the .general.bestassembly file to be the name and path of the sequence file
            metadata.general.bestassemblyfile = sample
            metadata.general.outputdirectory = os.path.join(self.outputdir, filename)
            metadata.general.logout = os.path.join(self.outputdir, filename, filename + '_out.txt')
            metadata.general.logerr = os.path.join(self.outputdir, filename, filename + '_err.txt')
            # Append the metadata for each sample to the list of samples
            self.runmetadata.samples.append(metadata)
