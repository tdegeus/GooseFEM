Data storage
============

Glossary
--------

+-----------+-------------------------------+-------------------------------------+
| Alias     | Shape                         | Description                         |
+===========+===============================+=====================================+
| *dofval*  | [ndof]                        | degrees of freedom                  |
+-----------+-------------------------------+-------------------------------------+
| *nodevec* | [nnode, ndim]                 | nodal vector                        |
+-----------+-------------------------------+-------------------------------------+
| *elemvec* | [nelem, nne, ndim]            | nodal vector stored per element     |
+-----------+-------------------------------+-------------------------------------+
| *elemmat* | [nelem, nne*ndim, nne*ndim]   | matrix stored per element           |
+-----------+-------------------------------+-------------------------------------+
| *qscalar* | [nelem, nip]                  | scalar stored per integration point |
+-----------+-------------------------------+-------------------------------------+
| *qtensor* | [nelem, nip, tdim, tdim, ...] | tensor stored per integration point |
+-----------+-------------------------------+-------------------------------------+

dofval
------

*   Degrees-of-freedom

*   Shape [ndof]

*   :code:`xt::xtensor<double,1>`

nodevec
-------

*   Nodal vectors

*   Shape [nnode, ndim]

*   :code:`xt::xtensor<double,2>`

elemvec
-------

*   Nodal vectors stored per element

*   Allows treatment of all elements independently, no connectivity needed

*   Shape [nelem, nne, ndim]

*   :code:`xt::xtensor<double,3>`

elemmat
-------

*   Matrices stored per element

*   Shape [nelem, nne*ndim, nne*ndim]

*   :code:`xt::xtensor<double,3>`

qscalar
-------

*   Scalars stored per integration point

*   Shape [nelem, nip]

*   :code:`xt::xtensor<double,2>`


qtensor (2nd order)
-------------------

*   2nd-order tensors stored per integration point

*   For certain elements, the number of dimensions of the tensor can be larger than the
    number of dimensions of the element (tdim >= ndim)

*   Shape [nelem, nip, tdim, tdim]

*   :code:`xt::xtensor<double,4>`

qtensor (4th order)
-------------------

*   4th-order tensors stored per integration point

*   For certain elements, the number of dimensions of the tensor can be larger than the
    number of dimensions of the element (tdim >= ndim)

*   Shape [nelem, nip, tdim, tdim, tdim, tdim]
*   :code:`xt::xtensor<double,4>`
