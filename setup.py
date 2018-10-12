from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

setup(name='drpclusters',
      version='1.1.1',
      description='Structural clustering of disulfide-rich peptides',
      url='https://github.com/dbarkan/drp-clusters/',
      author='Dave Barkan',
      author_email='davebarkan@gmail.com',
      license='GPL',
      #packages=['drpclusters'],
      packages=find_packages(),
      #include_package_data=True,
      #package_data={'drpclusters': ['exampleData/drp_list.txt']},
      zip_safe=False)
