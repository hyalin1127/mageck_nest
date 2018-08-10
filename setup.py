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
    version='3.0',
    author='Chen-Hao Chen, Wei Li',
    url='https://chenhaochen@bitbucket.org/liulab/mageck_nest.git',
    packages=['mageck_nest_v3'],
    package_dir={'mageck_nest_v3':'mageck_nest_v3'},
    package_data={'mageck_nest_v3':['*.p']}
  );


if __name__ == '__main__':
  main();
