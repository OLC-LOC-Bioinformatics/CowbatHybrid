#!/usr/bin/env python

import shutil
import logging
import subprocess


def check_dependencies():
    all_dependencies_good = True
    # Some stuff we don't care about version too much. For those, just check that the executable is present.
    dependencies = ['blastn',
                    'mob_recon',
                    'CLARK',
                    'NanoPlot',
                    'wtdbg2',
                    'wtpoa-cns',
                    'prodigal',
                    'sistr_cmd',
                    'mash',  # Not sure if screen functionality needed - if yes, update to require mash >=2.0
                    'GeneSeekr']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            logging.error('ERROR: Could not find dependency {}. Check that it is accessible from your $PATH'.format(dependency))
            all_dependencies_good = False
    # Other things have very specific versions - for those, need to actually check specific version.
    # wtdbg2 needs v2.1, 2.2 will require a change to the call in assembly.py
    # wtdbg2 doesn't get a version flag until 2.2. Not currently in bioconda - skip for now?

    return all_dependencies_good
