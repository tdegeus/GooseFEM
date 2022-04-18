****************
Getting GooseFEM
****************

Using conda
===========

The easiest is to use *conda* to install *GooseFEM*:

.. tabs::

    .. group-tab:: C++

        .. code-block:: bash

            conda install -c conda-forge goosefem

    .. group-tab:: Python

        .. code-block:: bash

            conda install -c conda-forge python-goosefem

        .. warning::

            This package does not benefit from (optional!) *xsimd* optimisation,
            as it is not compiled on your hardware.
            You'll have to compile by hand to benefit from *xsimd* optimisation.


This will install all the necessary runtime dependencies as well.

.. tip::

    The runtime dependencies (for both the C++ and the Python APIs)
    are also listed in ``environment.yaml``.
    One could install those dependencies in an activated environment by::

        conda env update --file environment.yaml

    This will install the dependencies to compile and run the code, tests, and examples.

From source
===========

.. tip::

    It could be instructive to see how configuration / compilation is done in the
    continuous integration: :download:`.github/workflows/ci.yml <../.github/workflows/ci.yml>`

.. tabs::

    .. group-tab:: C++

        Start by installing the dependencies, for example using *conda*::

            conda install -c conda-forge cmake xtensor

        Then, download the package::

            git checkout https://github.com/tdegeus/GooseFEM.git
            cd GooseFEM

        Install headers, *CMake* and *pkg-config* support:

        .. code-block:: none

            cmake -Bbuild
            cd build
            cmake --install .

        .. tip::

            To install in a loaded conda environment use

            .. code-block:: none

                cmake -Bbuild -DCMAKE_INSTALL_PREFIX:PATH="${CONDA_PREFIX}"
                cd build
                cmake --install .

    .. group-tab:: Python

        Start by installing the dependencies, for example using *conda*::

            conda install -c conda-forge xtensor-python eigen xsimd

        Note that *xsimd* is optional, but recommended.

        Then, download the package::

            git checkout https://github.com/tdegeus/GooseFEM.git
            cd GooseFEM

        Install the package using::

            python -m pip install . -v

        To use hardware optimisation (using *xsimd*) use instead::

            SKBUILD_CONFIGURE_OPTIONS="-DUSE_SIMD=1" python -m pip install . -v

.. note::

    The version is determined from the latest git tag, and possible commits since that tag.
    Python's ``setuptools_scm`` is used to this end.
    In case that you are not working from a clone of the repository you have to set
    the version manually, **before** configuring with CMake:

    .. tabs::

        .. group-tab:: Unix

            export SETUPTOOLS_SCM_PRETEND_VERSION="1.2.3"

        .. group-tab:: Windows

            set SETUPTOOLS_SCM_PRETEND_VERSION="1.2.3"

Docs
====

.. tip::

    It could be instructive to see how configuration / compilation is done in the
    continuous integration:
    :download:`.github/workflows/gh-pages.yml <../.github/workflows/gh-pages.yml>`

There are two kinds of docs:

1.  The current docs, generated using sphinx:

    .. code-block:: none

        cd docs
        make html

    Then open ``_build/html/index.html`` in your browser.

2.  The doxygen docs of the C++ API:

    .. code-block:: none

        cmake -Bbuild -DBUILD_DOCS=1
        cd build
        make html

    Then open ``html/index.html`` in your browser.
