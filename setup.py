#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="cowbathybrid",
    version="0.1.1",
    packages=find_packages(),
    scripts=['cowbat-hybrid-assembly.py'],
    author="Mathu Malar",
    author_email="Mathu.Malar@inspection.gc.ca",
    url="https://github.com/OLC-LOC-Bioinformatics/CowbatHybrid.git",
    install_requires=['olctools',
                      'geneseekr',
                      'sipprverse',
                      'seaborn',
                      'pandas',
                      'numpy',
                      'pysam',
                      'cowbat',
                      'genomeqaml']
)
