****************
Getting GooseFEM
****************

Using conda
===========

The easiest is to use *conda* to install *GooseFEM*::

    conda install -c conda-forge goosefem

This will install all the necessary runtime dependencies as well.

.. tip::

    The runtime dependencies (for both the C++ and the Python APIs) are also listed in
    ``environment.yaml``.
    One could install those dependencies in an activated environment by:

    .. code-block:: cpp

        conda env update --file environment.yaml

    In addition, one could further extend one's environment
    to also run the tests and the examples using:

    .. code-block:: cpp

        conda env update --file environment_test.yaml
        conda env update --file environment_examples.yaml

    Note that ``environment_test.yaml`` and ``environment_examples.yaml`` extend the environment.
    In each case one **also** has to install the dependencies in ``environment.yaml``.


From source
===========

Download the package::

    git checkout https://github.com/tdegeus/GooseFEM.git
    cd GooseFEM

Install headers, *CMake* and *pkg-config* support::

    cmake .
    make install

.. _install_python:

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

.. _install_python_source:

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
