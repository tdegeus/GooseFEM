
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

.. todo::

  Description, compile, run, and view instructions

:download:`main.cpp <examples/Statics/FixedDisplacements/LinearElasticity/main.cpp>`
:download:`CMakeLists.txt <examples/Statics/FixedDisplacements/LinearElasticity/CMakeLists.txt>`
:download:`plot.py <examples/Statics/FixedDisplacements/LinearElasticity/plot.py>`

:download:`main.cpp <examples/Statics/FixedDisplacements/LinearElasticity/main_manualPartition.cpp>`

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
