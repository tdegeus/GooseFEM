
*******
Statics
*******

.. tip::

  A compact reader covering the basic theory is available `here <https://github.com/tdegeus/GooseFEM/docs/theory/readme.pdf>`_.

Basic example
=============

Consider a uniform linear elastic bar that is extended by a uniform fixed displacement on both sides. This problem can be modelled and discretised using symmetry as show below. In this example we will furthermore assume that the bar is sufficiently thick in the out-of-plane direction that it can be modelled using two-dimensional plane strain.

.. image:: examples/Statics/FixedDisplacements/LinearElasticity/problem-sketch.svg
  :width: 300px
  :align: center

|

Include library
---------------

The first step is to include the header-only library:

.. literalinclude:: examples/Statics/FixedDisplacements/LinearElasticity/main.cpp
   :language: cpp
   :lines: 1-4
   :emphasize-lines: 1-2

Note that for this example we also make use of a material model (`GMatElastic <https://www.github.com/tdegeus/GMatElastic>`_) and a library to write (and read) HDF5 files (`HighFive <https://www.github.com/BlueBrain/HighFive>`_).

Define mesh
-----------

We will define a mesh using GooseFEM and extract the relevant information from it:

.. literalinclude:: examples/Statics/FixedDisplacements/LinearElasticity/main.cpp
   :language: cpp
   :lines: 11-28
   :emphasize-lines: 2

As observed the ``mesh`` is a class, that has methods to extract the relevant information such as the nodal coordinates (``coor``), the connectivity (``conn``), the degrees-of-freedom per node (``dofs``) and several node-sets that will be used to impose the sketched boundary conditions (``nodesLeft``, ``nodesRight``, ``nodesTop``, ``nodesBottom``).

Note that:

* The connectivity (``conn``) pro


and the degrees-of-freedom per node (``dofs``) as they connect the three data-representations in GooseFEM.



Define partitioning
-------------------

In particular, we will reorder such that

.. math::

  \texttt{u} =
  \begin{bmatrix}
    \texttt{u}_u \\
    \texttt{u}_p
  \end{bmatrix}

where the subscript :math:`u` and :math:`p` respectively denote *Unknown* and *Prescribed* degrees-of-freedom.


For our example, we will therefore first collect all degrees-of-freedom that are prescribed:

.. literalinclude:: examples/Statics/FixedDisplacements/LinearElasticity/main.cpp
   :language: cpp
   :lines: 33-38




.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <examples/Statics/FixedDisplacements/LinearElasticity/main.cpp>`
:download:`CMakeLists.txt <examples/Statics/FixedDisplacements/LinearElasticity/CMakeLists.txt>`
:download:`plot.py <examples/Statics/FixedDisplacements/LinearElasticity/plot.py>`

:download:`main.cpp <examples/Statics/FixedDisplacements/LinearElasticity/main_manualPartition.cpp>`

:download:`main.py <examples/Statics/FixedDisplacements/LinearElasticity/main.py>`

Partially periodic boundary conditions
======================================

.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <./examples/Statics/PeriodicPartial/LinearElasticity/main.cpp>`
:download:`CMakeLists.txt <./examples/Statics/PeriodicPartial/LinearElasticity/CMakeLists.txt>`
:download:`plot.py <./examples/Statics/PeriodicPartial/LinearElasticity/plot.py>`

Periodic boundary conditions
============================

.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <./examples/Statics/Periodic/LinearElasticity/main.cpp>`
:download:`CMakeLists.txt <./examples/Statics/Periodic/LinearElasticity/CMakeLists.txt>`
:download:`plot.py <./examples/Statics/Periodic/LinearElasticity/plot.py>`

Non-linear material behaviour
=============================

.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <./examples/Statics/Periodic/NonLinearElasticity/main.cpp>`
:download:`CMakeLists.txt <./examples/Statics/Periodic/NonLinearElasticity/CMakeLists.txt>`
:download:`plot.py <./examples/Statics/Periodic/NonLinearElasticity/plot.py>`

Non-linear & history dependent material behaviour
=================================================

.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <./examples/Statics/Periodic/ElastoPlasticity/main.cpp>`
:download:`CMakeLists.txt <./examples/Statics/Periodic/ElastoPlasticity/CMakeLists.txt>`
:download:`plot.py <./examples/Statics/Periodic/ElastoPlasticity/plot.py>`

Finite strain
=============

.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <./examples/Statics/Periodic/ElastoPlasticFiniteStrainSimo/main.cpp>`
:download:`CMakeLists.txt <./examples/Statics/Periodic/ElastoPlasticFiniteStrainSimo/CMakeLists.txt>`
:download:`plot.py <./examples/Statics/Periodic/ElastoPlasticFiniteStrainSimo/plot.py>`
