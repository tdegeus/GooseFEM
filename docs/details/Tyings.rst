.. _Tyings:

******
Tyings
******

| :download:`GooseFEM/TyingsPeriodic.h <../../include/GooseFEM/TyingsPeriodic.h>`
| :download:`GooseFEM/TyingsPeriodic.hpp <../../include/GooseFEM/TyingsPeriodic.hpp>`

Tyings::Periodic
================

Periodic nodal tyings: partition the system (renumber the DOFs) in the following order [iii, iid]: first the independent DOFs and the the dependent DOFs. If, in addition, independent DOFs are prescribed the partitioning is [iiu, iip, iid], where iii = [iiu, iip]: first the unknown and then the prescribed DOFs.

Tyings::Periodic::nnd()
-----------------------

Return the dependent DOF-numbers.

Tyings::Periodic::nni()
-----------------------

Return the independent DOF-numbers.

Tyings::Periodic::nnu()
-----------------------

Return the unknown DOF-numbers.

Tyings::Periodic::nnp()
-----------------------

Return the prescribed DOF-numbers.

Tyings::Periodic::dofs()
------------------------

Renumbered DOFs per node [nnode, ndim].

Tyings::Periodic::control()
---------------------------

Control DOF, that should be fixed [ndim].

Tyings::Periodic::iid()
-----------------------

Return the dependent DOFs.

Tyings::Periodic::iii()
-----------------------

Return the independent DOFs.

Tyings::Periodic::iiu()
-----------------------

Return the unknown DOFs.

Tyings::Periodic::iip()
-----------------------

Return the prescribed DOFs.

Tyings::Periodic::Cdi()
-----------------------

Return the tying matrix, such that

.. math::

  u_d = C_{di} u_i

In addition, the tying matrix in terms of the partitioned system can be obtained:

.. math::

  u_d = [C_{du}, C_{dp}]^T [u_u, u_p] = C_{du} u_u + C_{dp} u_p

Tyings::Periodic::Cdu()
-----------------------

Return the tying matrix, such that:

.. math::

  u_d = [C_{du}, C_{dp}]^T [u_u, u_p] = C_{du} u_u + C_{dp} u_p

Tyings::Periodic::Cdp()
-----------------------

Return the tying matrix, such that:

.. math::

  u_d = [C_{du}, C_{dp}]^T [u_u, u_p] = C_{du} u_u + C_{dp} u_p

Tyings::Control
===============

Add virtual control nodes to the system.

Tyings::Control::coor()
-----------------------

Nodal coordinates, including the virtual control nodes [nnode, ndim].

Tyings::Control::dofs()
-----------------------

DOFs, including the virtual control nodes [nnode, ndim].

Tyings::Control::controlDofs()
------------------------------

Virtual control DOFs [ndim, ndim].

Tyings::Control::controlNodes()
-------------------------------

Virtual control nodes [ndim].
