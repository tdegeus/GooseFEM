.. _Vector:

******
Vector
******

| :download:`GooseFEM/Vector.h <../../include/GooseFEM/Vector.h>`
| :download:`GooseFEM/Vector.hpp <../../include/GooseFEM/Vector.hpp>`
| :download:`GooseFEM/VectorPartitioned.h <../../include/GooseFEM/VectorPartitioned.h>`
| :download:`GooseFEM/VectorPartitioned.hpp <../../include/GooseFEM/VectorPartitioned.hpp>`
| :download:`GooseFEM/VectorPartitionedTyings.h <../../include/GooseFEM/VectorPartitionedTyings.h>`
| :download:`GooseFEM/VectorPartitionedTyings.hpp <../../include/GooseFEM/VectorPartitionedTyings.hpp>`

Vector
======

Vector definition allowing transforming between "dofval", "nodevec", and "elemvec" representations. See :ref:`conventions_vector`.

Vector::nelem()
---------------

Return the number of elements.

Vector::nne()
-------------

Return the number of nodes-per-element.

Vector::nnode()
---------------

Return the number of nodes.

Vector::ndim()
--------------

Return the number of dimensions.

Vector::ndof()
--------------

Return the number of DOFs.

Vector::dofs()
--------------

Return the DOF-numbers per node [nnode, ndim].

Vector::copy(...)
-----------------

Copy "nodevec" to "nodevec".

Vector::asDofs(...)
-------------------

Convert "nodevec" or "elemvec" to "dofval".

.. warning::

  Verify that you don't need "assembleDofs(...)"

Vector::asNode(...)
-------------------

Convert "dofval" or "elemvec" to "nodevec".

.. warning::

  Verify that you don't need "assembleNode(...)"

Vector::asElement(...)
----------------------

Convert "dofval" or "nodevec" to "elemvec".

Vector::assembleDofs(...)
-------------------------

Convert "nodevec" or "elemvec" to "dofval".

.. warning::

  Verify that you don't need "asDofs(...)"

Vector::assembleNode(...)
-------------------------

Convert "dofval" or "elemvec" to "nodevec".

.. warning::

  Verify that you don't need "asNode(...)"

Vector::AllocateDofval(...)
---------------------------

Allocate (and initialize) "dofval".

Vector::AllocateNodevec(...)
----------------------------

Allocate (and initialize) "nodevec".

Vector::AllocateElemvec(...)
----------------------------

Allocate (and initialize) "elemvec".

Vector::AllocateElemmat(...)
----------------------------

Allocate (and initialize) "elemmat".

VectorPartitioned
=================

Partitioned vector definition allowing transforming between "dofval", "nodevec", and "elemvec" representations. See :ref:`conventions_vector`. The partitioning is such that the DOFs are ordered as "[iiu, iip]" with "iiu" the unknown DOFs and "iip" the prescribed DOFs.

VectorPartitioned::nelem()
--------------------------

Return the number of elements.

VectorPartitioned::nne()
------------------------

Return the number of nodes-per-element.

VectorPartitioned::nnode()
--------------------------

Return the number of nodes.

VectorPartitioned::ndim()
-------------------------

Return the number of dimensions.

VectorPartitioned::ndof()
-------------------------

Return the number of DOFs.

VectorPartitioned::nnu()
------------------------

Return the number of unknown DOFs.

VectorPartitioned::nnp()
------------------------

Return the number of prescribed DOFs.

VectorPartitioned::dofs()
-------------------------

Return the DOF-numbers per node [nnode, ndim].

VectorPartitioned::iiu()
------------------------

Return the unknown DOF-numbers per node [nnu].

VectorPartitioned::iip()
------------------------

Return the prescribed DOF-numbers per node [nnp].

VectorPartitioned::copy(...)
----------------------------

Copy "nodevec" to "nodevec".

VectorPartitioned::copy_u(...)
------------------------------

Copy the unknown DOFs from a "nodevec" to the unknown DOFs from another "nodevec".

VectorPartitioned::copy_p(...)
------------------------------

Copy the prescribed DOFs from a "nodevec" to the prescribed DOFs from another "nodevec".

VectorPartitioned::asDofs(...)
------------------------------

Convert "nodevec" or "elemvec" to "dofval".

.. warning::

  Verify that you don't need "assembleDofs(...)"

VectorPartitioned::asDofs_u(...)
--------------------------------

Convert "nodevec" or "elemvec" to "dofval" and extract the unknown DOFs "iiu".

.. warning::

  Verify that you don't need "assembleDofs(...)"

VectorPartitioned::asDofs_p(...)
--------------------------------

Convert "nodevec" or "elemvec" to "dofval" and extract the prescribed DOFs "iip".

.. warning::

  Verify that you don't need "assembleDofs(...)"

VectorPartitioned::asNode(...)
------------------------------

Convert "dofval" or "elemvec" to "nodevec".

.. warning::

  Verify that you don't need "assembleNode(...)"

VectorPartitioned::asElement(...)
---------------------------------

Convert "dofval" or "nodevec" to "elemvec".

VectorPartitioned::assembleDofs(...)
------------------------------------

Convert "nodevec" or "elemvec" to "dofval".

.. warning::

  Verify that you don't need "asDofs(...)"

VectorPartitioned::assembleDofs_u(...)
--------------------------------------

Convert "nodevec" or "elemvec" to "dofval"  and extract the unknown DOFs "iiu".

.. warning::

  Verify that you don't need "asDofs(...)"

VectorPartitioned::assembleDofs_p(...)
--------------------------------------

Convert "nodevec" or "elemvec" to "dofval"  and extract the prescribed DOFs "iip".

.. warning::

  Verify that you don't need "asDofs(...)"

VectorPartitioned::assembleNode(...)
------------------------------------

Convert "dofval" or "elemvec" to "nodevec".

.. warning::

  Verify that you don't need "asNode(...)"

VectorPartitioned::AllocateDofval(...)
--------------------------------------

Allocate (and initialize) "dofval".

VectorPartitioned::AllocateNodevec(...)
---------------------------------------

Allocate (and initialize) "nodevec".

VectorPartitioned::AllocateElemvec(...)
---------------------------------------

Allocate (and initialize) "elemvec".

VectorPartitioned::AllocateElemmat(...)
---------------------------------------

Allocate (and initialize) "elemmat".

VectorPartitionedTyings
=======================

Partitioned vector definition with nodal tyings allowing transforming between "dofval", "nodevec", and "elemvec" representations. See :ref:`conventions_vector`. The partitioning is such that the DOFs are ordered as "[iiu, iip, iid]" with "iiu" the unknown DOFs and "iip" the prescribed DOFs and "iid" the dependent DOFs.

VectorPartitionedTyings::nelem()
--------------------------------

Return the number of elements.

VectorPartitionedTyings::nne()
------------------------------

Return the number of nodes-per-element.

VectorPartitionedTyings::nnode()
--------------------------------

Return the number of nodes.

VectorPartitionedTyings::ndim()
-------------------------------

Return the number of dimensions.

VectorPartitionedTyings::ndof()
-------------------------------

Return the number of DOFs.

VectorPartitionedTyings::nnu()
------------------------------

Return the number of unknown DOFs.

VectorPartitionedTyings::nnp()
------------------------------

Return the number of prescribed DOFs.

VectorPartitionedTyings::nni()
------------------------------

Return the number of independent DOFs.

VectorPartitionedTyings::nnd()
------------------------------

Return the number of dependent DOFs.

VectorPartitionedTyings::dofs()
-------------------------------

Return the DOF-numbers per node [nnode, ndim].

VectorPartitionedTyings::iiu()
------------------------------

Return the unknown DOF-numbers per node [nnu].

VectorPartitionedTyings::iip()
------------------------------

Return the prescribed DOF-numbers per node [nnp].

VectorPartitionedTyings::iii()
------------------------------

Return the independent DOF-numbers per node [nni].

VectorPartitionedTyings::iid()
------------------------------

Return the dependent DOF-numbers per node [nnd].

VectorPartitionedTyings::copy(...)
----------------------------------

Copy "nodevec" to "nodevec".

VectorPartitionedTyings::copy_u(...)
------------------------------------

Copy the unknown DOFs from a "nodevec" to the unknown DOFs from another "nodevec".

VectorPartitionedTyings::copy_p(...)
------------------------------------

Copy the prescribed DOFs from a "nodevec" to the prescribed DOFs from another "nodevec".

VectorPartitionedTyings::asDofs(...)
------------------------------------

Convert "nodevec" or "elemvec" to "dofval".

.. warning::

  Verify that you don't need "assembleDofs(...)"

VectorPartitionedTyings::asDofs_i(...)
--------------------------------------

Convert "nodevec" or "elemvec" to "dofval" and extract the independent DOFs "iii". Choose to apply the tyings:

.. math::

  u_i = C_{di}^T u_d

.. warning::

  Verify that you don't need "assembleDofs(...)"

VectorPartitionedTyings::asNode(...)
------------------------------------

Convert "dofval" or "elemvec" to "nodevec".

.. warning::

  Verify that you don't need "assembleNode(...)"

VectorPartitionedTyings::asElement(...)
---------------------------------------

Convert "dofval" or "nodevec" to "elemvec".

VectorPartitionedTyings::assembleDofs(...)
------------------------------------------

Convert "nodevec" or "elemvec" to "dofval".

.. warning::

  Verify that you don't need "asDofs(...)"

VectorPartitionedTyings::assembleNode(...)
------------------------------------------

Convert "dofval" or "elemvec" to "nodevec".

.. warning::

  Verify that you don't need "asNode(...)"

VectorPartitionedTyings::AllocateDofval(...)
--------------------------------------------

Allocate (and initialize) "dofval".

VectorPartitionedTyings::AllocateNodevec(...)
---------------------------------------------

Allocate (and initialize) "nodevec".

VectorPartitionedTyings::AllocateElemvec(...)
---------------------------------------------

Allocate (and initialize) "elemvec".

VectorPartitionedTyings::AllocateElemmat(...)
---------------------------------------------

Allocate (and initialize) "elemmat".
