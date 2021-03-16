****************
Getting GooseFEM
****************

Using conda
===========

The easiest is to use *conda* to install *GooseFEM*::

    conda install -c conda-forge goosefem

This will install all the necessary runtime dependencies as well.

.. tip::

    The runtime dependencies (for both the C++ and the Python APIs)
    are also listed in ``environment.yaml``.
    One could install those dependencies in an activated environment by:

    .. code-block:: cpp

        conda env update --file environment.yaml

    This will install the dependencies to run tests and examples.


From source
===========

.. tip::

    It could be instructive to see how configuration / compilation is done in the
    continuous integration: :download:`.github/workflows/ci.yml <../.github/workflows/ci.yml>`

Download the package::

    git checkout https://github.com/tdegeus/GooseFEM.git
    cd GooseFEM

Install headers, *CMake* and *pkg-config* support:

.. code-block:: none

    cmake -Bbuild
    cd build
    cmake --install .

.. note::

    The version is determined from the latest git tag, and possible commits since that tag.
    Internally Python's ``setuptools_scm`` is used to this end.
    In case that you are not working from a clone of the repository you have to set
    the version manually, **before** configuring with CMake:

    .. code-block:: none

        export SETUPTOOLS_SCM_PRETEND_VERSION="1.2.3"
        cmake -Bbuild
        cd build
        cmake --install .

    In Windows replace the first line with

    .. code-block:: none

        set SETUPTOOLS_SCM_PRETEND_VERSION="1.2.3"

.. tip::

    To install in a loaded conda environment use

    .. code-block:: none

        cmake -Bbuild -DCMAKE_INSTALL_PREFIX:PATH="${CONDA_PREFIX}"
        cd build
        cmake --install .


.. _install_python:

Python interface
================

.. tip::

    It could be instructive to see how configuration / compilation is done in the
    continuous integration: :download:`.github/workflows/ci.yml <../.github/workflows/ci.yml>`

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

.. _install_docs:

Docs
====

.. tip::

    It could be instructive to see how configuration / compilation is done in the
    continuous integration:
    :download:`.github/workflows/gh-pages.yml <../.github/workflows/gh-pages.yml>`

To build the docs there are two steps to be made:

1.  Extract the code documentation using doxygen:

    .. code-block:: none

        cmake -Bbuild -DBUILD_DOCS=1
        cd build
        make docs

2.  Build the docs using sphinx:

    .. code-block:: none

        cd docs
        make html
