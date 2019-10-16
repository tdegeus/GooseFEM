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

.. note::

  A solver has to be chosen, see :ref:`linear_solver`.

Matrix::nelem()
---------------

Return the number of elements.

Matrix::nne()
-------------

Return the number of nodes-per-element.

Matrix::nnode()
---------------

Return the number of nodes.

Matrix::ndim()
--------------

Return the number of dimensions.

Matrix::ndof()
--------------

Return the number of DOFs.

Matrix::dofs()
--------------

Return the DOF-numbers per node [nnode, ndim].

Matrix::assemble(...)
---------------------

Assemble matrix from element matrices stored as "elemmat".

Matrix::solve(...)
------------------

Solve linear system.

MatrixPartitioned
=================

Partitioned matrix definition.

.. note::

  A solver has to be chosen, see :ref:`linear_solver`.

MatrixPartitioned::nelem()
--------------------------

Return the number of elements.

MatrixPartitioned::nne()
------------------------

Return the number of nodes-per-element.

MatrixPartitioned::nnode()
--------------------------

Return the number of nodes.

MatrixPartitioned::ndim()
-------------------------

Return the number of dimensions.

MatrixPartitioned::ndof()
-------------------------

Return the number of DOFs.

MatrixPartitioned::nnu()
------------------------

Return the number of unknown DOFs.

MatrixPartitioned::nnp()
------------------------

Return the number of prescribed DOFs.

MatrixPartitioned::dofs()
-------------------------

Return the DOF-numbers per node [nnode, ndim].

MatrixPartitioned::iiu()
------------------------

Return the unknown DOF-numbers per node [nnu].

MatrixPartitioned::iip()
------------------------

Return the prescribed DOF-numbers per node [nnp].

MatrixPartitioned::assemble(...)
--------------------------------

Assemble matrix from element matrices stored as "elemmat".

MatrixPartitioned::solve(...)
-----------------------------

Solve linear system.

MatrixPartitioned::solve_u(...)
-------------------------------

Solve linear system.

MatrixPartitioned::reaction(...)
--------------------------------

Compute reaction forces (part of "b" that corresponds to "x_p").

MatrixPartitioned::reaction_p(...)
----------------------------------

Compute reaction forces (part of "b" that corresponds to "x_p").

MatrixPartitionedTyings
=======================

Partitioned matrix definition with nodal tyings.

.. note::

  A solver has to be chosen, see :ref:`linear_solver`.

MatrixPartitionedTyings::nelem()
--------------------------------

Return the number of elements.

MatrixPartitionedTyings::nne()
------------------------------

Return the number of nodes-per-element.

MatrixPartitionedTyings::nnode()
--------------------------------

Return the number of nodes.

MatrixPartitionedTyings::ndim()
-------------------------------

Return the number of dimensions.

MatrixPartitionedTyings::ndof()
-------------------------------

Return the number of DOFs.

MatrixPartitionedTyings::nnu()
------------------------------

Return the number of unknown DOFs.

MatrixPartitionedTyings::nnp()
------------------------------

Return the number of prescribed DOFs.

MatrixPartitionedTyings::nni()
------------------------------

Return the number of independent DOFs.

MatrixPartitionedTyings::nnd()
------------------------------

Return the number of dependent DOFs.

MatrixPartitionedTyings::dofs()
-------------------------------

Return the DOF-numbers per node [nnode, ndim].

MatrixPartitionedTyings::iiu()
------------------------------

Return the unknown DOF-numbers per node [nnu].

MatrixPartitionedTyings::iip()
------------------------------

Return the prescribed DOF-numbers per node [nnp].

MatrixPartitionedTyings::iii()
------------------------------

Return the independent DOF-numbers per node [nni].

MatrixPartitionedTyings::iid()
------------------------------

Return the dependent DOF-numbers per node [nnd].

MatrixPartitionedTyings::assemble(...)
--------------------------------------

Assemble matrix from element matrices stored as "elemmat".

MatrixPartitionedTyings::solve(...)
-----------------------------------

Solve linear system.

MatrixPartitionedTyings::solve_u(...)
-------------------------------------

Solve linear system.

MatrixDiagonal
==============

Diagonal matrix definition.

MatrixDiagonal::nelem()
-----------------------

Return the number of elements.

MatrixDiagonal::nne()
---------------------

Return the number of nodes-per-element.

MatrixDiagonal::nnode()
-----------------------

Return the number of nodes.

MatrixDiagonal::ndim()
----------------------

Return the number of dimensions.

MatrixDiagonal::ndof()
----------------------

Return the number of DOFs.

MatrixDiagonal::dofs()
----------------------

Return the DOF-numbers per node [nnode, ndim].

MatrixDiagonal::assemble(...)
-----------------------------

Assemble matrix from element matrices stored as "elemmat".

MatrixDiagonal::dot(...)
------------------------

Dot-product:

.. math::

  b_i = A_{ij} x_j

MatrixDiagonal::solve(...)
--------------------------

Solve linear system.

MatrixDiagonal::AsDiagonal(...)
-------------------------------

Return matrix as diagonal matrix (column)

MatrixDiagonalPartitioned
=========================

Diagonal and partitioned matrix definition.

MatrixDiagonalPartitioned::nelem()
----------------------------------

Return the number of elements.

MatrixDiagonalPartitioned::nne()
--------------------------------

Return the number of nodes-per-element.

MatrixDiagonalPartitioned::nnode()
----------------------------------

Return the number of nodes.

MatrixDiagonalPartitioned::ndim()
---------------------------------

Return the number of dimensions.

MatrixDiagonalPartitioned::ndof()
---------------------------------

Return the number of DOFs.

MatrixDiagonalPartitioned::nnu()
--------------------------------

Return the number of unknown DOFs.

MatrixDiagonalPartitioned::nnp()
--------------------------------

Return the number of prescribed DOFs.

MatrixDiagonalPartitioned::dofs()
---------------------------------

Return the DOF-numbers per node [nnode, ndim].

MatrixDiagonalPartitioned::iiu()
--------------------------------

Return the unknown DOF-numbers per node [nnu].

MatrixDiagonalPartitioned::iip()
--------------------------------

Return the prescribed DOF-numbers per node [nnp].

MatrixDiagonalPartitioned::assemble(...)
----------------------------------------

Assemble matrix from element matrices stored as "elemmat".

MatrixDiagonalPartitioned::dot(...)
-----------------------------------

Dot-product:

.. math::

  b_i = A_{ij} x_j

MatrixDiagonalPartitioned::dot_u(...)
-------------------------------------

Dot-product:

.. math::

  b_i = A_{ij} x_j

MatrixDiagonalPartitioned::dot_p(...)
-------------------------------------

Dot-product:

.. math::

  b_i = A_{ij} x_j

MatrixDiagonalPartitioned::solve(...)
-------------------------------------

Solve linear system.

MatrixDiagonalPartitioned::solve_u(...)
---------------------------------------

Solve linear system.

MatrixDiagonalPartitioned::reaction(...)
----------------------------------------

Compute reaction forces (part of "b" that corresponds to "x_p").

MatrixDiagonalPartitioned::reaction_p(...)
------------------------------------------

Compute reaction forces (part of "b" that corresponds to "x_p").

MatrixDiagonalPartitioned::AsDiagonal(...)
------------------------------------------

Return matrix as diagonal matrix (column)

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

