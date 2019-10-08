
.. _compile:

*********
Compiling
*********

Introduction
============

This module is header only. So one just has to ``#include <GooseFEM/GooseFEM.h>``. somewhere in the source code, and to tell the compiler where the header-files are. For the latter, several ways are described below.

.. tip::

  Optimisation of crucial importance if you do not want to wait forever. Please `use the strategies provided by xtensor <https://xtensor.readthedocs.io/en/latest/build-options.html>`_. In particular, it is highly advice to use `xsimd <https://github.com/QuantStack/xsimd>`_ in addition to the usual optimisation flags.

.. note::

  This code depends on `eigen3 <https://github.com/RLovelett/eigen>`_ and `xtensor <https://github.com/QuantStack/xtensor>`_. Both are also header-only libraries. Both can be 'installed' using identical steps as described below.

Manual compiler flags
=====================

GNU / Clang
-----------

Add the following compiler's arguments (in addition to the arguments to include `xtensor <https://github.com/QuantStack/xtensor>`_):

.. code-block:: bash

  -I${PATH_TO_GOOSEFEM}/src

.. note:: **(Not recommended)**

  If you want to avoid separately including the header files using a compiler flag, ``git submodule`` is a nice way to go:

  1.  Include this module as a submodule using ``git submodule add https://github.com/tdegeus/GooseFEM.git``.

  2.  Replace the first line of this example by ``#include "GooseFEM/include/GooseFEM/GooseFEM.h"``.

      *If you decide to manually copy the header file, you might need to modify this relative path to your liking.*

  Or see :ref:`compile_automatic`. You can also combine the ``git submodule`` with any of the below compiling strategies.

.. _compile_automatic:

(Semi-)Automatic compiler flags
===============================

Install
-------

To enable (semi-)automatic build, one should 'install' ``GooseFEM`` somewhere.

Install system-wide (root)
^^^^^^^^^^^^^^^^^^^^^^^^^^

1.  Proceed to a (temporary) build directory. For example

    .. code-block:: bash

      $ cd /path/to/GooseFEM/src/build

2.  'Build' ``GooseFEM``

    .. code-block:: bash

      $ cmake ..
      $ make install

    (If you've used another build directory, change the first command to ``$ cmake /path/to/GooseFEM/src``)

Install in custom location (user)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1.  Proceed to a (temporary) build directory. For example

    .. code-block:: bash

      $ cd /path/to/GooseFEM/src/build

2.  'Build' ``GooseFEM`` to install it in a custom location

    .. code-block:: bash

      $ mkdir /custom/install/path
      $ cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/custom/install/path
      $ make install

    (If you've used another build directory, change the first command to ``$ cmake /path/to/GooseFEM/src``)

3.  Add the following path to your ``~/.bashrc`` (or ``~/.zshrc``):

    .. code-block:: bash

      export PKG_CONFIG_PATH=/custom/install/path/share/pkgconfig:$PKG_CONFIG_PATH

.. note:: **(Not recommended)**

  If you do not wish to use ``CMake`` for the installation, or you want to do something custom. You can of course. Follow these steps:

  1.  Copy the file ``include/GooseFEM.pc.in`` to ``GooseFEM.pc`` to some location that can be found by ``pkg_config`` (for example by adding ``export PKG_CONFIG_PATH=/path/to/GooseFEM.pc:$PKG_CONFIG_PATH`` to the ``.bashrc``).

  2.  Modify the line ``prefix=@CMAKE_INSTALL_PREFIX@`` to ``prefix=/path/to/GooseFEM``.

  3.  Modify the line ``Cflags: -I${prefix}/@INCLUDE_INSTALL_DIR@`` to ``Cflags: -I${prefix}/src``.

  4.  Modify the line ``Version: @GOOSEFEM_VERSION_NUMBER@`` to reflect the correct release version.

Compiler arguments from 'pkg-config'
------------------------------------

Instead of ``-I...`` one can now use

.. code-block:: bash

  `pkg-config --cflags GooseFEM` -std=c++14

as compiler argument.

Compiler arguments from 'cmake'
-------------------------------

Add the following to your ``CMakeLists.txt``:

.. code-block:: cmake

  set(CMAKE_CXX_STANDARD 14)

  find_package(PkgConfig)

  pkg_check_modules(GOOSEFEM REQUIRED GooseFEM)
  include_directories(${GOOSEFEM_INCLUDE_DIRS})
