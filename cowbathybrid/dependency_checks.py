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
                    'prodigal',
                    'pilon',
                    'porechop',
                    'sistr',
                    'mash',  # Not sure if screen functionality needed - if yes, update to require mash >=2.0
                    'GeneSeekr',
                    'famap',
                    'fahash',
                    'classify.py']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            logging.error('ERROR: Could not find dependency {}. Check that it is accessible from your $PATH'.format(dependency))
            all_dependencies_good = False
        else:
            logging.debug('Found {} at {}'.format(dependency, shutil.which(dependency)))
    # Other things have very specific versions - for those, need to actually check specific version.
    # Unicycler version 0.4.7 stalls seemingly at random. As far as I can tell 0.4.4 does not suffer from the same
    # issue, so we'll enforce 0.4.4
    try:
        unicycler_version = subprocess.check_output('unicycler --version', shell=True).decode('utf-8').split()[1]
    except subprocess.CalledProcessError:
        unicycler_version = 'Not Found'
    if unicycler_version != 'v0.4.4':
        logging.error('ERROR: Unicycler version found was {}, but this pipeline requires v0.4.4 - please install '
                      'the correct version and try again.'.format(unicycler_version))
        all_dependencies_good = False
    return all_dependencies_good
