
********************
Notes for developers
********************

Testing
=======

Please add relevant tests to ``test/``. To run all tests:

.. code-block:: cpp

    mkdir build/
    cmake .. -DBUILD_TESTS=ON
    ./test/test

Building the docs
=================

To set-up build dependencies consider using Conda

.. code-block:: cpp

    conda env update --file environment_docs.yaml

Then,

.. code-block:: cpp

    cd docs
    doxygen
    make html

.. note::

    By default the latest version of the Python API on conda-forge is used,
    this may results in outdated docs.
    If you want the latests docs, please build and install the Python module first,
    see :ref:`install_python_source`.
