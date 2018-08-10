#!/usr/bin/env python
'''
MAGeCK set up script
'''
from __future__ import print_function;

import os
import sys
from distutils.core import setup, Extension
from subprocess import call as subpcall
from distutils.command.install import install as DistutilsInstall

try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

def main():
  # check python version
  if float(sys.version[:3])<3.0:
    sys.stderr.write("CRITICAL: Python version must be >=3.0!\n")
    sys.exit(1);

  setup(
    name='mageck_nest',
    author='Chen-Hao Chen, Wei Li',
    author_email='hyalin1127@gmail.com',
    url='https://hyalin1127@bitbucket.org/hyalin1127/mageck_nest.git',
    packages=['mageck_nest'],
    scripts = ['bin/mageck_nest'],
    package_dir={'mageck_nest':'mageck_nest'},
    #package_data={'mageck_nest':['*.p']}
  );

if __name__ == '__main__':
  main();
