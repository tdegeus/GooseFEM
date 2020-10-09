****************
Getting GooseFEM
****************


Using conda
===========

The easiest is to use *conda* to install *GooseFEM*::

    conda install -c conda-forge goosefem


From source
===========

Download the package::

    git checkout https://github.com/tdegeus/GooseFEM.git
    cd GooseFEM

Install headers, *CMake* and *pkg-config* support::

    cmake .
    make install


Python interface
================

Using conda
^^^^^^^^^^^

The quickest (but not the most efficient!) is to use *conda* to install *GooseFEM*::

    conda install -c conda-forge python-goosefem

.. warning::

    This package does not benefit from *xsimd* optimisation,
    as it is not compiled on your hardware.
    You'll have to compile by hand to benefit from *xsimd* optimisation.

From source
^^^^^^^^^^^

Start by installing the dependencies, for example using *conda*::

    conda install -c conda-forge pyxtensor eigen xsimd

Note that *xsimd* is optional, but recommended.

.. note::

    You can also use::

        python -m pip install pyxtensor pybind11

    for use without *conda*. Note that you install *Eigen* and *xsimd* yourself
    in such a way that Python can find it in order to use it.

Then, download the package::

    git checkout https://github.com/tdegeus/GooseFEM.git
    cd GooseFEM

Install the package using::

    python -m pip install .

.. note::

    The following will give more readable output::

        python setup.py build
        python setup.py install
