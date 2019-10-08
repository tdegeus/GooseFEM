
**************
Linear elastic
**************

Consider a uniform linear elastic bar that is extended by a uniform fixed displacement on both sides. This problem can be modelled and discretised using symmetry as show below. In this example we will furthermore assume that the bar is sufficiently thick in the out-of-plane direction that it can be modelled using two-dimensional plane strain.

.. image:: statics/FixedDisplacements_LinearElastic/problem-sketch.svg
  :width: 300px
  :align: center

|

Include library
===============

The first step is to include the header-only library:

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 1-4
   :emphasize-lines: 1-2

Note that for this example we also make use of a material model (`GMatElastic <https://www.github.com/tdegeus/GMatElastic>`_) and a library to write (and read) HDF5 files (`HighFive <https://www.github.com/BlueBrain/HighFive>`_).

Define mesh
===========

We will define a mesh using GooseFEM and extract the relevant information from it:

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 11-28
   :emphasize-lines: 2

As observed the "mesh" is a class, that has methods to extract the relevant information such as the nodal coordinates ("coor"), the connectivity ("conn"), the degrees-of-freedom per node ("dofs") and several node-sets that will be used to impose the sketched boundary conditions ("nodesLeft", "nodesRight", "nodesTop", "nodesBottom").

Note that:

* The connectivity ("conn") contains information of which nodes, in which order, belong to which element.
* The degrees-of-freedom per node ("dofs") contains information of how a nodal vector (a vector stored per node) can be transformed to a list of degrees-of-freedom as used in the linear system (although this can be mostly done automatically as we will see below).

.. seealso::

  :ref:`conventions_terminology`

Define partitioning
===================

In particular, we will reorder such that degrees-of-freedom are ordered such that

.. math::

  \texttt{u} =
  \begin{bmatrix}
    \texttt{u}_u \\
    \texttt{u}_p
  \end{bmatrix}

where the subscript :math:`u` and :math:`p` respectively denote *Unknown* and *Prescribed* degrees-of-freedom.

For our example, we will therefore first collect all degrees-of-freedom that are prescribed:

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 33-38




.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <statics/FixedDisplacements_LinearElastic/example/main.cpp>`
:download:`CMakeLists.txt <statics/FixedDisplacements_LinearElastic/example/CMakeLists.txt>`
:download:`plot.py <statics/FixedDisplacements_LinearElastic/example/plot.py>`

:download:`main.cpp <statics/FixedDisplacements_LinearElastic/manual_partition//main.cpp>`

:download:`main.py <statics/FixedDisplacements_LinearElastic/example/main.py>`

:download:`main.py <statics/FixedDisplacements_LinearElastic/manual_partition//main.py>`
