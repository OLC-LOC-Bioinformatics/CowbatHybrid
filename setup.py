#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="cowbathybrid",
    version="0.1.2",
    packages=find_packages(),
    scripts=['cowbat-hybrid-assembly.py'],
    author="Mathu MalarAndrew Low",
    author_email=" Mathu.Malar@inspection.gc.ca, andrew.low@canada.ca",
    url="https://github.com/OLC-LOC-Bioinformatics/CowbatHybrid",
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
