.. _MeshHex8:

**********
Mesh::Hex8
**********

| :download:`GooseFEM/MeshHex8.h <../../include/GooseFEM/MeshHex8.h>`
| :download:`GooseFEM/MeshHex8.hpp <../../include/GooseFEM/MeshHex8.hpp>`

Naming convention
=================

The following naming convention is used:

* **Front**: all nodes whose coordinates :math:`0 \leq x \leq L_x`, :math:`0 \leq y \leq L_y`, :math:`z = 0`.
* **Back**: all nodes whose coordinates :math:`0 \leq x \leq L_x`, :math:`0 \leq y \leq L_y`, :math:`z = L_z`.
* **Bottom**: all nodes whose coordinates :math:`0 \leq x \leq L_x`, :math:`0 \leq z \leq L_z`, :math:`y = 0`.
* **Top**: all nodes whose coordinates :math:`0 \leq x \leq L_x`, :math:`0 \leq z \leq L_z`, :math:`y = L_y`.
* **Left**: all nodes whose coordinates :math:`0 \leq y \leq L_y`, :math:`0 \leq z \leq L_z`, :math:`x = 0`.
* **Right**: all nodes whose coordinates :math:`0 \leq y \leq L_y`, :math:`0 \leq z \leq L_z`, :math:`x = L_x`.

The edges and corners follow from the intersections, i.e.

* **FrontBottomEdge**: all nodes whose coordinates :math:`0 \leq x \leq L_x`, :math:`y = 0`, :math:`z = 0`.
* ...
* **FrontBottomLeftCorner**: the node whose coordinate :math:`x = 0`, :math:`y = 0`, :math:`z = 0`.
* ...

.. image:: figures/MeshHex8/naming_convention.svg
  :width: 350px
  :align: center

Mesh::Hex8::Regular
===================

Regular mesh.

.. todo::

  Describe, illustrate

.. todo::

  More methods and classes

.. todo::

  figures/MeshHex8/Regular/example...

.. todo::

  figures/MeshHex8/Regular/nodes...


Mesh::Hex8::FineLayer
=====================

Mesh with a middle plane that is fine the middle, and becomes course away from this plane.

.. todo::

  Describe, illustrate

.. todo::

  figures/MeshHex8/FineLayer/example...

.. todo::

  figures/MeshHex8/FineLayer/nodes...
