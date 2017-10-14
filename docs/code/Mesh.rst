
******************
<GooseFEM/Mesh.h>
******************

``GooseFEM::Mesh::dofs``
========================

Get a sequential list of DOF-numbers for each vector-component of each node.

.. code-block:: cpp

  MatS GooseFEM::Mesh::dofs ( size_t nnode , size_t ndim )

For example for 3 nodes in 2 dimensions the output is

.. math::

  \begin{bmatrix}
    0 & 1 \\
    2 & 3 \\
    4 & 5
  \end{bmatrix}

``GooseFEM::Mesh::renumber``
============================

* Renumber indices to lowest possible indices (returns a copy, input not modified).

  .. code-block:: cpp

    MatS GooseFEM::Mesh::renumber ( const MatS &in )

  For example:

  .. math::

    \begin{bmatrix}
      0 & 1 \\
      0 & 1 \\
      5 & 4
    \end{bmatrix}

  is renumbered to

  .. math::

    \begin{bmatrix}
      0 & 1 \\
      0 & 1 \\
      3 & 2
    \end{bmatrix}

* Reorder indices such that some items are at the beginning or the end (returns a copy, input not modified).

  .. code-block:: cpp

    MatS GooseFEM::Mesh::renumber ( const MatS &in , const VecS &idx, std::string location="end" );

  For example:

  .. math::

    \mathrm{in} =
    \begin{bmatrix}
      0 & 1 \\
      0 & 1 \\
      3 & 2 \\
      4 & 5
    \end{bmatrix}

  with

  .. math::

    \mathrm{idx} =
    \begin{bmatrix}
      6 & 4
    \end{bmatrix}

  Implies that ``in`` is renumbered such that the 6th item becomes the one-before-last item (:math:`5 \rightarrow 4`), and the 4th item become the last (:math:`3 \rightarrow 5`). The remaining items are renumbered to the lowest index while keeping the same order. The result:

  .. math::

    \begin{bmatrix}
      0 & 1 \\
      0 & 1 \\
      5 & 2 \\
      4 & 3
    \end{bmatrix}






