
********************
GooseFEM::Mesh::Hex8
********************

[:download:`MeshHex8.h <../src/GooseFEM/MeshHex8.h>`, :download:`MeshHex8.cpp <../src/GooseFEM/MeshHex8.cpp>`]

Naming convention
=================

.. todo::

  Update

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

``GooseFEM::Mesh::Hex8::Regular``
=================================

.. todo::

  Update


``GooseFEM::Mesh::Hex8::FineLayer``
===================================

.. todo::

  Update


.. .. code-block:: cpp

..   GooseFEM::Mesh::Quad4::Regular(size_t nelx, size_t nely, double h=1.);

.. Regular mesh of linear quadrilaterals in two-dimensions. The element edges are all of the same size :math:`h` (by default equal to one), optional scaling can be applied afterwards. For example the mesh shown below that consists of 21 x 11 elements. In that image the element numbers are indicated with a color, and likewise for the boundary nodes.



.. Methods:

.. .. code-block:: cpp

..   // A matrix with on each row a nodal coordinate:
..   // [ x , y ]
..   MatD = GooseFEM::Mesh::Quad4::Regular.coor();

..   // A matrix with the connectivity, with on each row to the nodes of each element
..   MatS = GooseFEM::Mesh::Quad4::Regular.conn();

..   // A list of boundary nodes
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesBottom();
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesTop();
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesLeft();
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesRight();

..   // A matrix with periodic node pairs on each row:
..   // [ independent nodes, dependent nodes ]
..   MatS = GooseFEM::Mesh::Quad4::Regular.nodesPeriodic();

..   // The node at the origin
..   size_t = GooseFEM::Mesh::Quad4::Regular.nodeOrigin();

..   // A matrix with DOF-numbers: two per node in sequential order
..   MatS = GooseFEM::Mesh::Quad4::Regular.dofs();

..   // A matrix with DOF-numbers: two per node in sequential order
..   // All the periodic repetitions are eliminated from the system
..   MatS = GooseFEM::Mesh::Quad4::Regular.dofsPeriodic();

.. ``GooseFEM::Mesh::Quad4::FineLayer``
.. ====================================

.. Regular mesh with a fine layer of quadrilateral elements, and coarser elements above and below.

.. .. image:: figures/MeshQuad4/FineLayer/example.svg
..   :width: 500px
..   :align: center

.. .. note::

..   The coarsening depends strongly on the desired number of elements in horizontal elements. The becomes clear from the following example:

..   .. code-block:: cpp

..     mesh = GooseFEM::Mesh::Quad4::FineLayer(6*9  ,51); // left   image :  546 elements
..     mesh = GooseFEM::Mesh::Quad4::FineLayer(6*9+3,51); // middle image :  703 elements
..     mesh = GooseFEM::Mesh::Quad4::FineLayer(6*9+1,51); // right  image : 2915 elements

..   .. image:: figures/MeshQuad4/FineLayer/behavior.svg
..     :width: 1000px
..     :align: center

.. Methods:

.. .. code-block:: cpp

..   // A matrix with on each row a nodal coordinate:
..   // [ x , y ]
..   MatD = GooseFEM::Mesh::Quad4::Regular.coor();

..   // A matrix with the connectivity, with on each row to the nodes of each element
..   MatS = GooseFEM::Mesh::Quad4::Regular.conn();

..   // A list of boundary nodes
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesBottom();
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesTop();
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesLeft();
..   ColS = GooseFEM::Mesh::Quad4::Regular.nodesRight();

..   // A matrix with periodic node pairs on each row:
..   // [ independent nodes, dependent nodes ]
..   MatS = GooseFEM::Mesh::Quad4::Regular.nodesPeriodic();

..   // The node at the origin
..   size_t = GooseFEM::Mesh::Quad4::Regular.nodeOrigin();

..   // A matrix with DOF-numbers: two per node in sequential order
..   MatS = GooseFEM::Mesh::Quad4::Regular.dofs();

..   // A matrix with DOF-numbers: two per node in sequential order
..   // All the periodic repetitions are eliminated from the system
..   MatS = GooseFEM::Mesh::Quad4::Regular.dofsPeriodic();

..   // A list with the element numbers of the fine elements in the center of the mesh
..   // (highlighted in the plot below)
..   ColS = GooseFEM::Mesh::Quad4::FineLayer.elementsFine();

..     .. image:: figures/MeshQuad4/FineLayer/example_elementsFine.svg
..       :width: 500px
..       :align: center
