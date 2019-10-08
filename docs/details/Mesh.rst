
****
Mesh
****

| :download:`GooseFEM/Mesh.h <../../include/GooseFEM/Mesh.h>`
| :download:`GooseFEM/Mesh.hpp <../../include/GooseFEM/Mesh.hpp>`

Mesh::dofs
----------

Get a sequential list of DOF-numbers for each vector-component of each node. For example for 3 nodes in 2 dimensions the output is

.. math::

  \begin{bmatrix}
    0 & 1 \\
    2 & 3 \\
    4 & 5
  \end{bmatrix}

Mesh::Renumber
--------------

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

  One can use the wrapper function "GooseFEM::reorder" or the class "Mesh::Reorder" to get more advanced features.

Mesh::Reorder
-------------

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

  One can use the wrapper function "GooseFEM::reorder" or the class "Mesh::Reorder" to get more advanced features.

Mesh::coordination
------------------

Get the number of elements connected to each node.

Mesh::elem2node
---------------

Get the element numbers (columns) that are connected to each node (rows).
