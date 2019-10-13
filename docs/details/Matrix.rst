.. _Matrix:

******
Matrix
******

| :download:`GooseFEM/Matrix.h <../../include/GooseFEM/Matrix.h>`
| :download:`GooseFEM/Matrix.hpp <../../include/GooseFEM/Matrix.hpp>`
| :download:`GooseFEM/MatrixPartitioned.h <../../include/GooseFEM/MatrixPartitioned.h>`
| :download:`GooseFEM/MatrixPartitioned.hpp <../../include/GooseFEM/MatrixPartitioned.hpp>`
| :download:`GooseFEM/MatrixPartitionedTyings.h <../../include/GooseFEM/MatrixPartitionedTyings.h>`
| :download:`GooseFEM/MatrixPartitionedTyings.hpp <../../include/GooseFEM/MatrixPartitionedTyings.hpp>`
| :download:`GooseFEM/MatrixDiagonal.h <../../include/GooseFEM/MatrixDiagonal.h>`
| :download:`GooseFEM/MatrixDiagonal.hpp <../../include/GooseFEM/MatrixDiagonal.hpp>`
| :download:`GooseFEM/MatrixDiagonalPartitioned.h <../../include/GooseFEM/MatrixDiagonalPartitioned.h>`
| :download:`GooseFEM/MatrixDiagonalPartitioned.hpp <../../include/GooseFEM/MatrixDiagonalPartitioned.hpp>`

Matrix
======

Matrix definition.

.. todo::

  Complete.

MatrixPartitioned
=================

Partitioned matrix definition.

.. todo::

  Complete.

MatrixPartitionedTyings
=======================

Partitioned matrix definition with nodal tyings.

.. todo::

  Complete.

MatrixDiagonal
==============

Diagonal matrix definition.

.. todo::

  Complete.

MatrixDiagonalPartitioned
=========================

Diagonal and partitioned matrix definition.

.. todo::

  Complete.

.. _linear_solver:

Linear solver
=============

The classes ``GooseFEM:::MatrixPartitioned`` and ``GooseFEM:::MatrixPartitionedTyings`` make use of a library to solver the linear system (stored as a sparse matrix). In particular the Eigen library and its plug-ins are used. To use the library's default solver:

.. code-block:: cpp

    #include <Eigen/Eigen>
    #include <GooseFEM/GooseFEM.h>

    int main()
    {
        ...

        GooseFEM::MatrixPartitioned<> K(...);

        ...

        return 0;
    }

The default solver can, however, be quite slow. Therefore Eigen has quite some `plug-ins <http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html>`_ for the solver. GooseFEM allows the use of Eigen's Sparse Solver Concept to use such plug-ins. For example, to use SuiteSparse:

.. code-block:: cpp

    #include <Eigen/Eigen>
    #include <Eigen/CholmodSupport>
    #include <GooseFEM/GooseFEM.h>

    int main()
    {
        ...

        GooseFEM::MatrixPartitioned<Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>> K(...);

        ...

        return 0;
    }

.. todo::

    1.  `Download SuiteSparse <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_.

    2.  Extract the downloaded ``SuiteSparse-X.X.X.tar.gz```.

    3.  Compile the library by:

        .. code-block:: bash

            cd /path/to/SuiteSparse
            make

    .. code-block:: bash

        clang++ -I/path/to/include/eigen3 -I/path/to/lapack/include -L/path/to/lapack/lib -I/path/to/openblas/include -L/path/to/openblas/lib -lopenblas -I/path/to/SuiteSparse/include -L/path/to/SuiteSparse/lib -lumfpack -lamd -lcholmod -lsuitesparseconfig -lm -std=c++14 -Wall -Wextra -pedantic -march=native -O3  -o example example.cpp

