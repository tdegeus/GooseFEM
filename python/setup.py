desc = '''
GooseFEM is a C++ module, wrapped in Python, that provides several predefined finite element meshes.
The original C++ module also includes element definitions and several standard finite element
simulations.
'''

from setuptools import setup, Extension

import sys,re
import setuptools
import pybind11
import cppmat

header = open('../src/GooseFEM/GooseFEM.h','r').read()
world  = re.split('(.*)(\#define GOOSEFEM_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split('(.*)(\#define GOOSEFEM_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split('(.*)(\#define GOOSEFEM_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'GooseFEM',
    ['python_interface.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      cppmat  .get_include(False),
      cppmat  .get_include(True ),
      cppmat  .find_eigen()
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
  install_requires = ['pybind11>=2.2.0','cppmat>=0.2.15'],
  cmdclass         = {'build_ext': cppmat.BuildExt},
  zip_safe         = False,
)
