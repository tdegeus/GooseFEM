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

Mesh::Hex8::Regular::nelem()
----------------------------

Return number of elements.

Mesh::Hex8::Regular::nnode()
----------------------------

Return number of nodes.

Mesh::Hex8::Regular::nne()
--------------------------

Return number of nodes-per-element (= 3).

Mesh::Hex8::Regular::ndim()
---------------------------

Return number of dimensions (= 2).

Mesh::Hex8::Regular::getElementType()
-------------------------------------

Return element-type.

Mesh::Hex8::Regular::coor()
---------------------------

Return nodal coordinates [nnode, ndim].

Mesh::Hex8::Regular::conn()
---------------------------

Return connectivity [nelem, nne].

Mesh::Hex8::Regular::nodesXXXEdge()
-----------------------------------

Node numbers along the "Bottom", "Top", "Left", or "Right" edge.

Mesh::Hex8::Regular::nodesXXXOpenEdge()
---------------------------------------

Node numbers along the "Bottom", "Top", "Left", or "Right" edge, excluding the corners.

Mesh::Hex8::Regular::nodesXXXCorner()
-------------------------------------

Node number of one of the corners (e.g. "BottomLeft").

Mesh::Hex8::Regular::nodesPeriodic()
------------------------------------

Periodic node pairs. Each row contains on pair of (independent, dependent) node numbers. The output shape is thus [n_pairs, 2].

Mesh::Hex8::Regular::nodesOrigin()
----------------------------------

Bottom-left node, used as reference for periodicity.

Mesh::Hex8::Regular::dofs()
---------------------------

DOF-numbers for each component of each node (sequential). The output shape is thus [nnode, ndim].

Mesh::Hex8::Regular::dofsPeriodic()
-----------------------------------

DOF-numbers for each component of each node, for the case that the periodicity if fully eliminated. The output shape is thus [nnode, ndim].

Mesh::Hex8::Regular::elementgrid()
----------------------------------

Return element numbers as matrix [nelz, nely, nelx].

Mesh::Hex8::FineLayer
=====================

Mesh with a fine layer in the middle, and that becomes course away from this plane (see image below). Note coursening can only be done if the number of elements in x- and y-direction is dividable by 3, and that it is only optimal if the number of elements in x- and y-direction is a factor of 3. Note that the number of elements in the vertical direction is specified as the number of times the unit element (the number of times "h" the height should be), and that this number is only a target: the algorithm chooses in accordance with the applied coursing.

Mesh::Hex8::FineLayer::nelem()
------------------------------

Return number of elements.

Mesh::Hex8::FineLayer::nnode()
------------------------------

Return number of nodes.

Mesh::Hex8::FineLayer::nne()
----------------------------

Return number of nodes-per-element (= 3).

Mesh::Hex8::FineLayer::ndim()
-----------------------------

Return number of dimensions (= 2).

Mesh::Hex8::FineLayer::nelx()
-----------------------------

Number of elements in horizontal direction (along the weak layer) (matches input).

Mesh::Hex8::FineLayer::nely()
-----------------------------

Actual number of elements unit elements in vertical direction (actual number of times "h" the mesh is heigh).

Mesh::Hex8::FineLayer::h()
--------------------------

Unit edge size (matches input).

Mesh::Hex8::FineLayer::getElementType()
---------------------------------------

Return element-type.

Mesh::Hex8::FineLayer::coor()
-----------------------------

Return nodal coordinates [nnode, ndim].

Mesh::Hex8::FineLayer::conn()
-----------------------------

Return connectivity [nelem, nne].

Mesh::Hex8::FineLayer::nodesXXXEdge()
-------------------------------------

Node numbers along the "Bottom", "Top", "Left", "Right", "Front", or "Back" edge.

Mesh::Hex8::FineLayer::nodesXXXOpenEdge()
-----------------------------------------

Node numbers along the "Bottom", "Top", "Left", "Right", "Front", or "Back" edge, excluding the corners.

Mesh::Hex8::FineLayer::nodesXXXCorner()
---------------------------------------

Node number of one of the corners (e.g. "FrontBottomLeft").

Mesh::Hex8::FineLayer::nodesPeriodic()
--------------------------------------

Periodic node pairs. Each row contains on pair of (independent, dependent) node numbers. The output shape is thus [n_pairs, 2].

Mesh::Hex8::FineLayer::nodesOrigin()
------------------------------------

Bottom-left node, used as reference for periodicity.

Mesh::Hex8::FineLayer::dofs()
-----------------------------

DOF-numbers for each component of each node (sequential). The output shape is thus [nnode, ndim].

Mesh::Hex8::FineLayer::dofsPeriodic()
-------------------------------------

DOF-numbers for each component of each node, for the case that the periodicity if fully eliminated. The output shape is thus [nnode, ndim].

Mesh::Hex8::FineLayer::elementsMiddleLayer()
--------------------------------------------

Element numbers of the middle, fine, layer
