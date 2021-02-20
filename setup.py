
from setuptools import setup, Extension

import re
import os
import pybind11
import pyxtensor
from os import environ

version = environ.get('PKG_VERSION')

if version is None:
    from setuptools_scm import get_version
    version = get_version()

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

build.c_opts['unix'] += ['-DGOOSEFEM_VERSION="{0:s}"'.format(version)]
build.c_opts['msvc'] += ['/DGOOSEFEM_VERSION=\\"{0:s}\\"'.format(version)]

ext_modules = [Extension(
    'GooseFEM',
    ['python/main.cpp'],
    include_dirs = include_dirs,
    language = 'c++')]

setup(
    name = 'GooseFEM',
    description = 'Finite element meshes, quadrature, and assembly tools',
    long_description = 'Finite element meshes, quadrature, and assembly tools',
    version = version,
    license = 'GPLv3',
    author = 'Tom de Geus',
    author_email = 'tom@geus.me',
    url = 'https://github.com/tdegeus/GooseFEM',
    ext_modules = ext_modules,
    install_requires = ['pybind11', 'pyxtensor'],
    cmdclass = {'build_ext': build},
    zip_safe = False)
