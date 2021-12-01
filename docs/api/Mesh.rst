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

.. doxygenfunction:: GooseFEM::Mesh::edgesize(const C&, const E&, ElementType)
    :project: GooseFEM

edgesize
--------

.. doxygenfunction:: GooseFEM::Mesh::edgesize(const C&, const E&)
    :project: GooseFEM

centers
-------

.. doxygenfunction:: GooseFEM::Mesh::centers(const C&, const E&, ElementType)
    :project: GooseFEM

.. doxygenfunction:: GooseFEM::Mesh::centers(const C&, const E&)
    :project: GooseFEM

elemmap2nodemap
---------------

.. doxygenfunction:: GooseFEM::Mesh::elemmap2nodemap(const T&, const C&, const E&, ElementType)
    :project: GooseFEM

.. doxygenfunction:: GooseFEM::Mesh::elemmap2nodemap(const T&, const C&, const E&)
    :project: GooseFEM

center_of_gravity
-----------------

.. doxygenfunction:: GooseFEM::Mesh::center_of_gravity(const C&, const E&, ElementType)
    :project: GooseFEM

.. doxygenfunction:: GooseFEM::Mesh::center_of_gravity(const C&, const E&)
    :project: GooseFEM
