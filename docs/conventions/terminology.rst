.. _conventions_terminology:

Terminology
===========

Glossary
--------

+---------+---------+-------------------------------------+
| Alias   | Type    | Description                         |
+=========+=========+=====================================+
| nnode   | size    | number of nodes                     |
+---------+---------+-------------------------------------+
| ndim    | size    | number of dimensions                |
+---------+---------+-------------------------------------+
| nelem   | size    | number of elements                  |
+---------+---------+-------------------------------------+
| nne     | size    | number of nodes per element         |
+---------+---------+-------------------------------------+
| tdim    | size    | number of dimensions of tensor      |
+---------+---------+-------------------------------------+
| ---     | ---     | ---                                 |
+---------+---------+-------------------------------------+
| dofval  | array   | degrees of freedom                  |
+---------+---------+-------------------------------------+
| nodevec | array   | nodal vector                        |
+---------+---------+-------------------------------------+
| elemvec | array   | nodal vector stored per element     |
+---------+---------+-------------------------------------+
| elemmat | array   | matrix stored per element           |
+---------+---------+-------------------------------------+
| qscalar | array   | scalar stored per integration point |
+---------+---------+-------------------------------------+
| qtensor | array   | tensor stored per integration point |
+---------+---------+-------------------------------------+
| ---     | ---     | ---                                 |
+---------+---------+-------------------------------------+
| coor    | name    | nodal coordinates                   |
+---------+---------+-------------------------------------+
| conn    | name    | connectivity                        |
+---------+---------+-------------------------------------+
| ---     | ---     | ---                                 |
+---------+---------+-------------------------------------+
| Tri3    | element | 2-d triangular element (3 nodes)    |
+---------+---------+-------------------------------------+
| Quad4   | element | 2-d quadrilateral element (4 nodes) |
+---------+---------+-------------------------------------+
| Hex8    | element | 3-d hexahedral element (8 nodes)    |
+---------+---------+-------------------------------------+

Coordinates
-----------

* Nodal coordinates: each node has one row
* Shape :math:`\left[ n_\text{node} \times n_\text{dim} \right]`
* Denoted: *coor*

Connectivity
------------

* Node numbers per element: each element has one row
* Shape :math:`\left[ n_\text{elem} \times n_\text{nodes-per-element} \right]`
* Denoted: *conn*

Sizes
-----

+-------+---------------------------------+------------------------------------+
| Alias | Description                     |                                    |
+=======+=================================+====================================+
| nnode | number of nodes                 | :math:`n_\text{nodes}`             |
+-------+---------------------------------+------------------------------------+
| ndim  | number of dimensions            | :math:`n_\text{ndim}`              |
+-------+---------------------------------+------------------------------------+
| nelem | number of elements              | :math:`n_\text{elements}`          |
+-------+---------------------------------+------------------------------------+
| nne   | number of nodes per element     | :math:`n_\text{nodes-per-element}` |
+-------+---------------------------------+------------------------------------+
| tdim  | number of dimensions of tensor  | :math:`n_\text{ndim}`              |
+-------+---------------------------------+------------------------------------+

