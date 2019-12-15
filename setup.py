desc = '''
GooseFEM is a C++ module, wrapped in Python, that provides several predefined finite element meshes.
The original C++ module also includes element definitions and several standard finite element
simulations.
'''

from setuptools import setup, Extension

import re
import os
import pybind11
import pyxtensor

header = open('include/GooseFEM/config.h','r').read()
major = re.split(r'(.*)(\#define GOOSEFEM_VERSION_MAJOR\ )([0-9]+)(.*)', header)[3]
minor = re.split(r'(.*)(\#define GOOSEFEM_VERSION_MINOR\ )([0-9]+)(.*)', header)[3]
patch = re.split(r'(.*)(\#define GOOSEFEM_VERSION_PATCH\ )([0-9]+)(.*)', header)[3]

__version__ = '.'.join([major, minor, patch])

include_dirs = [
    os.path.abspath('include/'),
    pyxtensor.find_pyxtensor(),
    pyxtensor.find_pybind11(),
    pyxtensor.find_xtensor(),
    pyxtensor.find_xtl(),
    pyxtensor.find_eigen()]

build = pyxtensor.BuildExt

xsimd = pyxtensor.find_xsimd()

if xsimd:
    if len(xsimd) > 0:
        include_dirs += [xsimd]
        build.c_opts['unix'] += ['-march=native', '-DXTENSOR_USE_XSIMD']
        build.c_opts['msvc'] += ['/DXTENSOR_USE_XSIMD']

ext_modules = [Extension(
    'GooseFEM',
    ['python/main.cpp'],
    include_dirs = include_dirs,
    language = 'c++')]

setup(
    name = 'GooseFEM',
    description = 'Finite element meshes, quadrature, and assembly tools',
    long_description = desc,
    version = __version__,
    license = 'GPLv3',
    author = 'Tom de Geus',
    author_email = 'tom@geus.me',
    url = 'https://github.com/tdegeus/GooseFEM',
    ext_modules = ext_modules,
    install_requires = ['pybind11>=2.2.0', 'pyxtensor>=0.1.1'],
    cmdclass = {'build_ext': build},
    zip_safe = False)
