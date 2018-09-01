desc = '''
GooseFEM is a C++ module, wrapped in Python, that provides several predefined finite element meshes.
The original C++ module also includes element definitions and several standard finite element
simulations.
'''

from setuptools import setup, Extension

import sys, re
import setuptools
import pybind11
import pyxtensor

header = open('include/xGooseFEM/GooseFEM.h','r').read()
world  = re.split(r'(.*)(\#define XGOOSEFEM_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split(r'(.*)(\#define XGOOSEFEM_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split(r'(.*)(\#define XGOOSEFEM_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'GooseFEM',
    ['include/xGooseFEM/python.cpp'],
    include_dirs=[
      pybind11 .get_include(False),
      pybind11 .get_include(True ),
      pyxtensor.get_include(False),
      pyxtensor.get_include(True ),
      pyxtensor.find_xtensor(),
      pyxtensor.find_xtl(),
      pyxtensor.find_eigen(),
    ],
    language='c++'
  ),
]

setup(
  name             = 'GooseFEM',
  description      = 'Finite element meshes, quadrature, and assembly tools',
  long_description = desc,
  version          = __version__,
  license          = 'GPLv3',
  author           = 'Tom de Geus',
  author_email     = 'tom@geus.me',
  url              = 'https://github.com/tdegeus/GooseFEM',
  ext_modules      = ext_modules,
  install_requires = ['pybind11>=2.2.0','pyxtensor>=0.0.1'],
  cmdclass         = {'build_ext': pyxtensor.BuildExt},
  zip_safe         = False,
)
