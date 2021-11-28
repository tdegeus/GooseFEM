.. _conventions_terminology:

Terminology
===========

Sizes
-----

+-----------+----------------------------------+
| Alias     | Description                      |
+===========+==================================+
| *nnode*   | number of nodes                  |
+-----------+----------------------------------+
| *ndim*    | number of dimensions             |
+-----------+----------------------------------+
| *nelem*   | number of elements               |
+-----------+----------------------------------+
| *nne*     | number of nodes per element      |
+-----------+----------------------------------+
| *tdim*    | number of dimensions of a tensor |
+-----------+----------------------------------+

Arrays
------

+-----------+-------------------------------------+
| Alias     | Description                         |
+===========+=====================================+
| *dofval*  | degrees of freedom                  |
+-----------+-------------------------------------+
| *nodevec* | nodal vector                        |
+-----------+-------------------------------------+
| *elemvec* | nodal vector stored per element     |
+-----------+-------------------------------------+
| *elemmat* | matrix stored per element           |
+-----------+-------------------------------------+
| *qscalar* | scalar stored per integration point |
+-----------+-------------------------------------+
| *qtensor* | tensor stored per integration point |
+-----------+-------------------------------------+

Names
-----

+-----------+-------------------+
| Alias     | Description       |
+===========+===================+
| *coor*    | nodal coordinates |
+-----------+-------------------+
| *conn*    | connectivity      |
+-----------+-------------------+

Elements
--------

+-----------+-------------------------------------+
| Alias     | Description                         |
+===========+=====================================+
| *Tri3*    | 2-d triangular element (3 nodes)    |
+-----------+-------------------------------------+
| *Quad4*   | 2-d quadrilateral element (4 nodes) |
+-----------+-------------------------------------+
| *Hex8*    | 3-d hexahedral element (8 nodes)    |
+-----------+-------------------------------------+

Coordinates
-----------

* Nodal coordinates: each node has one row
* Shape [*nnode*, *ndim*]
* Denoted: *coor*

Connectivity
------------

* Node numbers per element: each element has one row
* Shape [*nelem*, *nne*]
* Denoted: *conn*
