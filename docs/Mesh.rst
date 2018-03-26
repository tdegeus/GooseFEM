
**************
GooseFEM::Mesh
**************

[:download:`GooseFEM/Mesh.h <../src/GooseFEM/Mesh.h>`, :download:`GooseFEM/Mesh.cpp <../src/GooseFEM/Mesh.cpp>`]

``GooseFEM::Mesh::dofs``
========================

.. code-block:: cpp

  GooseFEM::MatS GooseFEM::Mesh::dofs(size_t nnode, size_t ndim)

Get a sequential list of DOF-numbers for each vector-component of each node. For example for 3 nodes in 2 dimensions the output is

.. math::

  \begin{bmatrix}
    0 & 1 \\
    2 & 3 \\
    4 & 5
  \end{bmatrix}

``GooseFEM::Mesh::renumber``
============================

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

.. note:: Generic interface

  .. code-block:: cpp

    template<class InputIterator, class OutputIterator>
    void renumber(const InputIterator first, const InputIterator last, const OutputIterator result)

``GooseFEM::Mesh::reorder``
===========================

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

.. note:: Generic interface

  .. code-block:: cpp

    template<class InputIterator, class OutputIterator, class IndexIterator>
    void reorder(const InputIterator first, const InputIterator last, const OutputIterator result, const IndexIterator first_index, const IndexIterator last_index, std::string location)

``GooseFEM::Mesh::elem2node``
=============================

.. code-block:: cpp

  GooseFEM::SpMatS GooseFEM::Mesh::elem2node(const GooseFEM::MatS &conn)

Return a sparse matrix which contains the element numbers (columns) that are connected to each node (rows).

.. warning::

  One should not confuse the element ``0`` when this matrix is converted to a dense matrix. When this is done all the 'missing' items are filled in as zero, which does have a meaning here.


