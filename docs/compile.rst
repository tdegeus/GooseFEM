
.. _compile:

***************
Compiling (C++)
***************

Introduction
============

This module is header only. So one just has to ``#include <GooseFEM/GooseFEM.h>``, and to tell the compiler where the header-files are.

Using CMake
===========

The following structure can be used for your project's ``CMakeLists.txt``:

.. code-block:: cmake

    find_package(GooseFEM REQUIRED)

    add_executable(myexec mymain.cpp)

    target_link_libraries(myexec PRIVATE
        GooseFEM
        xtensor::optimize
        xtensor::use_xsimd)

See the `documentation of xtensor <https://xtensor.readthedocs.io/en/latest/>`_ concerning optimisation.

.. note::

    There are additional targets available to expedite your ``CMakeLists.txt``:

    *   ``GooseFEM::assert``: enable GooseFEM assertions.
    *   ``GooseFEM::debug``: enable GooseFEM and xtensor assertions (slow).
    *   ``GooseFEM::compiler_warnings``: enable compiler warnings.

By hand
=======

Presuming that the compiler is ``c++``, compile using::

    c++ -I/path/to/GooseFEM/include ...

Note that you have to take care of the *xtensor* and *Eigen* dependencies, the C++ version, optimisation, enabling *xsimd*, ...

Using pkg-config
================

Find the location of the headers can be automatised using *pkg-config*::

    pkg-config --cflags GooseFEM
