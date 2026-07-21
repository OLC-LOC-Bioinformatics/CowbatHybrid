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
                    'flye',
                    'classify.py',
                    'filtlong']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            logging.error('ERROR: Could not find dependency {}. Check that it is accessible from your $PATH'.format(dependency))
            all_dependencies_good = False
        else:
            logging.debug('Found {} at {}'.format(dependency, shutil.which(dependency)))

    # Unicycler version check
    try:
        unicycler_version = subprocess.check_output(
            'unicycler --version', shell=True
        ).decode('utf-8').split()[1]
    except subprocess.CalledProcessError:
        unicycler_version = 'Not Found'

    if unicycler_version not in ('v0.4.4', 'v0.5.1'):
        logging.error(
            'ERROR: Unicycler version found was {}, but this pipeline requires one of: v0.4.4, v0.5.1 - '
            'please install a supported version and try again.'.format(unicycler_version)
        )
        all_dependencies_good = False

    return all_dependencies_good
