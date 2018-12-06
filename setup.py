#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="cowbathybrid",
    version="0.1.1",
    packages=find_packages(),
    scripts=['cowbathybrid/cowbat-hybrid-assembly.py'],
    author="Andrew Low",
    author_email="andrew.low@canada.ca",
    url="https://github.com/lowandrew/cowbat-hybrid",
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
