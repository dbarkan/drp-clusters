from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

setup(name='drpclusters',
      version='0.1.1',
      description='Structural clustering of disulfide-rich peptides',
      url='https://github.com/dbarkan/drp-clusters/',
      author='Dave Barkan',
      author_email='davebarkan@gmail.com',
      install_requires=['biopython>=1.7.2'],
      license='GPL',
      packages=['drpclusters'],
      #include_package_data=True,
      #package_data={'drpclusters': ['data/longerFraction.txt']},
      zip_safe=False)
