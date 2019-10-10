.. _Element:

*******
Element
*******

| :download:`GooseFEM/Element.h <../../include/GooseFEM/Element.h>`
| :download:`GooseFEM/Element.hpp <../../include/GooseFEM/Element.hpp>`

.. tip::

  To take the gradients and integral with respect to updated coordinates (i.e. to do updated Lagrange), use the ".update_x(...)" method to update the nodal coordinates and re-evaluate the shape function gradients and integration volumes.

Element::asElementVector
========================

Convert nodal vector "[nnode, ndim]" to nodal vector stored per element "[nelem, nne, ndim]".

.. todo::

  Describe.

Element::assembleNodeVector
===========================

Assemble nodal vector stored per element "[nelem, nne, ndim]" to nodal vector "[nnode, ndim]".

.. todo::

  Describe.

Element::isSequential
=====================

Check that DOFs leave no holes.

.. todo::

  Describe.

Element::isDiagonal
===================

Check structure of the matrices stored per element "[nelem, nne*ndim, nne*ndim]".

.. todo::

  Describe.
