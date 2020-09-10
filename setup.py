#!/usr/bin/env python3

from setuptools import setup, find_packages
from crispr_read_counts.version import version

config = {
  'name': 'crispr-read-counts',
  'description': 'Code to count reads for CRISPR',
  'author': 'Yaobo Xu',
  'url': 'https://gitlab.internal.sanger.ac.uk/CancerIT/crisprReadCountsPy',
  'download_url': '',
  'author_email': 'cgphelp@sanger.ac.uk',
  'version': version,
  'python_requires': '>= 3.6',
  'setup_requires': ['pytest'],
  'install_requires': [
    'click==7.1.2',
    'pysam'],
  'packages': find_packages(),
  'entry_points': {
    'console_scripts': ['crisprReadCounts=crispr_read_counts.command_line:main'],
  },
}

setup(**config)
