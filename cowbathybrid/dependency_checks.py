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
                    'unicycler',
                    'prodigal',
                    'sistr',
                    'mash',  # Not sure if screen functionality needed - if yes, update to require mash >=2.0
                    'GeneSeekr']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            logging.error('ERROR: Could not find dependency {}. Check that it is accessible from your $PATH'.format(dependency))
            all_dependencies_good = False
    # Other things have very specific versions - for those, need to actually check specific version.

    return all_dependencies_good
