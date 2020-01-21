****************
Getting GooseFEM
****************


Compiling
=========

Getting GooseFEM
----------------

Using conda
^^^^^^^^^^^

The easiest is to use *conda* to install *GooseFEM*::

    conda install -c conda-forge goosefem

From source
^^^^^^^^^^^

Download the package::

    git checkout https://github.com/tdegeus/GooseFEM.git
    cd GooseFEM

Install headers, *CMake* and *pkg-config* support::

    cmake .
    make install

Compiling
---------

Using CMake
^^^^^^^^^^^

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
^^^^^^^

Presuming that the compiler is ``c++``, compile using::

    c++ -I/path/to/GooseFEM/include ...

Note that you have to take care of the *xtensor* and *Eigen* dependencies, the C++ version, optimisation, enabling *xsimd*, ...

Using pkg-config
^^^^^^^^^^^^^^^^

Find the location of the headers can be automatised using *pkg-config*::

    pkg-config --cflags GooseFEM


Python interface
================

Using conda
^^^^^^^^^^^

The quickest (but not the most efficient!) is to use *conda* to install *GooseFEM*::

    conda install -c conda-forge python-goosefem

.. warning::

    This package does not benefit from *xsimd* optimisation, as it is not compiled on your hardware. Therefore compiling by hand is recommended.

From source
^^^^^^^^^^^

Start by installing the dependencies, for example using *conda*::

    conda install -c conda-forge pyxtensor eigen xsimd

Note that *xsimd* is optional, but recommended.

.. note::

    You can also use::

        python -m pip install pyxtensor pybind11

    for use without *conda*. Note that you install *Eigen* and *xsimd* yourself in such a way that Python can find it in order to use it.

Then, download the package::

    git checkout https://github.com/tdegeus/GooseFEM.git
    cd GooseFEM

Install the package using::

    python -m pip install .

.. note::

    The following will give more readable output::

        python setup.py build
        python setup.py install
