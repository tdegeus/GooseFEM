****
Mesh
****

Element-type
============

ElementType
-----------

.. doxygenenum:: GooseFEM::Mesh::ElementType
   :project: GooseFEM

defaultElementType
------------------

.. doxygenfunction:: GooseFEM::Mesh::defaultElementType
    :project: GooseFEM

Stitch meshes
=============

Stitch
------

.. doxygenclass:: GooseFEM::Mesh::Stitch
   :project: GooseFEM
   :members:

ManualStitch
------------

.. doxygenclass:: GooseFEM::Mesh::ManualStitch
   :project: GooseFEM
   :members:

overlapping
-----------

.. doxygenfunction:: GooseFEM::Mesh::overlapping
    :project: GooseFEM

Renumber lists
==============

Renumber
--------

.. doxygenclass:: GooseFEM::Mesh::Renumber
   :project: GooseFEM
   :members:

renumber
--------

.. doxygenfunction:: GooseFEM::Mesh::renumber
    :project: GooseFEM


Reorder
-------

.. doxygenclass:: GooseFEM::Mesh::Reorder
   :project: GooseFEM
   :members:

Basic mesh properties
=====================

dofs
----

.. doxygenfunction:: GooseFEM::Mesh::dofs
    :project: GooseFEM

coordination
------------

.. doxygenfunction:: GooseFEM::Mesh::coordination
    :project: GooseFEM

elem2node
---------

.. doxygenfunction:: GooseFEM::Mesh::elem2node
    :project: GooseFEM

edgesize
--------

.. doxygenfunction:: GooseFEM::Mesh::edgesize(const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&, ElementType)
    :project: GooseFEM

edgesize
--------

.. doxygenfunction:: GooseFEM::Mesh::edgesize(const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&)
    :project: GooseFEM

centers
-------

.. doxygenfunction:: GooseFEM::Mesh::centers(const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&, ElementType)
    :project: GooseFEM

.. doxygenfunction:: GooseFEM::Mesh::centers(const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&)
    :project: GooseFEM

elemmap2nodemap
---------------

.. doxygenfunction:: GooseFEM::Mesh::elemmap2nodemap(const xt::xtensor<size_t, 1>&, const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&, ElementType)
    :project: GooseFEM

.. doxygenfunction:: GooseFEM::Mesh::elemmap2nodemap(const xt::xtensor<size_t, 1>&, const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&)
    :project: GooseFEM
