.. _Element:

*******
Element
*******

| :download:`GooseFEM/Element.h <../../include/GooseFEM/Element.h>`
| :download:`GooseFEM/Element.hpp <../../include/GooseFEM/Element.hpp>`

Element::asElementVector
========================

Convert nodal vector "[nnode, ndim]" to nodal vector stored per element "[nelem, nne, ndim]".

Element::assembleNodeVector
===========================

Assemble nodal vector stored per element "[nelem, nne, ndim]" to nodal vector "[nnode, ndim]".

Element::isSequential
=====================

Check that DOFs leave no holes.

Element::isDiagonal
===================

Check structure of the matrices stored per element "[nelem, nne*ndim, nne*ndim]" to be diagonal (check that all off-diagonal entries have a value lower than a small numerical tolerance).
