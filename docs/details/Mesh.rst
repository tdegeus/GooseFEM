
****
Mesh
****

Generic methods
===============

[:download:`GooseFEM/Mesh.h <../../include/GooseFEM/Mesh.h>`, :download:`GooseFEM/Mesh.hpp <../../include/GooseFEM/Mesh.hpp>`]

GooseFEM::Mesh::dofs
--------------------

.. code-block:: cpp

  GooseFEM::MatS GooseFEM::Mesh::dofs(size_t nnode, size_t ndim)

Get a sequential list of DOF-numbers for each vector-component of each node. For example for 3 nodes in 2 dimensions the output is

.. math::

  \begin{bmatrix}
    0 & 1 \\
    2 & 3 \\
    4 & 5
  \end{bmatrix}

GooseFEM::Mesh::renumber
------------------------

.. code-block:: cpp

  GooseFEM::MatS GooseFEM::Mesh::renumber(const GooseFEM::MatS &dofs)

Renumber (DOF) indices to lowest possible indices. For example:

.. math::

  \begin{bmatrix}
    0 & 1 \\
    5 & 4
  \end{bmatrix}

is renumbered to

.. math::

  \begin{bmatrix}
    0 & 1 \\
    3 & 2
  \end{bmatrix}

Or, in pseudo-code, the result of this function is that:

.. code-block:: python

  dofs = renumber(dofs)

  sort(unique(dofs[:])) == range(max(dofs+1))

.. tip::

  A generic interface using iterator is available if you do not which to use the default Eigen interface.

GooseFEM::Mesh::reorder
-----------------------

.. code-block:: cpp

  GooseFEM::MatS GooseFEM::Mesh::reorder(const GooseFEM::MatS &dofs, const ColS &idx, std::string location="end")

Reorder (DOF) indices such to the lowest possible indices, such that some items are at the beginning or the end. For example:

.. math::

  \mathrm{dofs} =
  \begin{bmatrix}
    0 & 1 \\
    2 & 3 \\
    4 & 5
  \end{bmatrix}

with

.. math::

  \mathrm{idx} =
  \begin{bmatrix}
    0 & 1
  \end{bmatrix}

Implies that ``dofs`` is renumbered such that 0 becomes the one-before-last index (:math:`0 \rightarrow 4`), and the 1 becomes the last index (:math:`1 \rightarrow 5`). The remaining items are renumbered to the lowest index while keeping the same order. The result:

.. math::

  \begin{bmatrix}
    4 & 5 \\
    0 & 1 \\
    2 & 3
  \end{bmatrix}

.. tip::

  A generic interface using iterator is available if you do not which to use the default Eigen interface.

GooseFEM::Mesh::elem2node
-------------------------

.. code-block:: cpp

  GooseFEM::SpMatS GooseFEM::Mesh::elem2node(const GooseFEM::MatS &conn)

Return a sparse matrix which contains the element numbers (columns) that are connected to each node (rows).

.. warning::

  One should not confuse the element ``0`` when this matrix is converted to a dense matrix. When this is done all the 'missing' items are filled in as zero, which does have a meaning here.

Predefined meshes
=================

GooseFEM::Mesh::Tri3
--------------------

[:download:`GooseFEM/MeshTri3.h <../../include/GooseFEM/MeshTri3.h>`, :download:`GooseFEM/MeshTri3.hpp <../../include/GooseFEM/MeshTri3.hpp>`]

GooseFEM::Mesh::Tri3::Regular
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.


GooseFEM::Mesh::Hex8
--------------------

[:download:`MeshHex8.h <../../include/GooseFEM/MeshHex8.h>`, :download:`MeshHex8.hpp <../../include/GooseFEM/MeshHex8.hpp>`]

Naming convention
^^^^^^^^^^^^^^^^^

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

GooseFEM::Mesh::Hex8::Regular
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regular mesh.

GooseFEM::Mesh::Hex8::FineLayer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mesh with a middle plane that is fine the middle, and becomes course away from this plane.

Type specific methods
=====================

GooseFEM::Mesh::Tri3
--------------------

GooseFEM::Mesh::Tri3::Regular
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[:download:`GooseFEM/MeshTri3.h <../../include/GooseFEM/MeshTri3.h>`, :download:`GooseFEM/MeshTri3.hpp <../../include/GooseFEM/MeshTri3.hpp>`]


GooseFEM::Mesh::Tri3::getOrientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.

GooseFEM::Mesh::Tri3::setOrientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.

GooseFEM::Mesh::Tri3::retriangulate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.

GooseFEM::Mesh::Tri3::TriUpdate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.

GooseFEM::Mesh::Tri3::Edge
^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.
