
********
Dynamics
********

Inertial dynamics, no velocity dependence
=========================================

.. todo::

  Description and example

Verlet
------

Basic structure
^^^^^^^^^^^^^^^

.. code-block:: cpp

  // position, velocity, acceleration (and history: last increment)
  xt::xtensor<double,2> u   = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v   = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> a   = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v_n = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> a_n = xt::zeros<double>({nnode, ndim});

  // residual force
  xt::xtensor<double,2> fr = xt::zeros<double>({nnode, ndim});

  // compute mass matrix
  // (often assumed constant & diagonal, remove either assumption if needed)
  GooseFEM::MatrixDiagonal M(...);

  ...

  // time increments
  for ( ... )
  {
    // store history
    xt::noalias(v_n) = v;
    xt::noalias(a_n) = a;

    // new displacement
    xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt,2.) * a;

    // new residual force (and mass matrix if needed)
    ...

    // new acceleration
    M.solve(fr, a);

    // new velocity
    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );
  }

Full example
^^^^^^^^^^^^

.. todo::

  Compile, run, and view instructions

:download:`main.cpp <./examples/Dynamics/ElasticLaminate-Verlet/main.cpp>`
:download:`CMakeLists.txt <./examples/Dynamics/ElasticLaminate-Verlet/CMakeLists.txt>`
:download:`plot.py <./examples/Dynamics/ElasticLaminate-Verlet/plot.py>`

.. code-block:: bash

  mkdir build
  cd build
  make
  ./main
  python3 plot.py

Inertial dynamics, velocity dependence
======================================

.. todo::

  Description and example

Velocity Verlet
---------------

.. code-block:: cpp

  // position, velocity, acceleration (and history: last increment)
  xt::xtensor<double,2> u   = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v   = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> a   = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v_n = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> a_n = xt::zeros<double>({nnode, ndim});

  // residual force
  xt::xtensor<double,2> fr = xt::zeros<double>({nnode, ndim});

  // compute mass matrix
  // (often assumed constant & diagonal, remove either assumption if needed)
  GooseFEM::MatrixDiagonal M(...);

  ...

  // time increments
  for ( ... )
  {
    // store history
    xt::noalias(v_n) = v;
    xt::noalias(a_n) = a;

    // new displacement
    xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt,2.) * a;

    // update residual force (and mass matrix if needed)
    ...

    // estimate 1: new velocity
    xt::noalias(v) = v_n + dt * a_n;

    // estimate 1: new residual force (and mass matrix if needed)
    ...

    // estimate 1: new acceleration
    M.solve(fr, a);

    // estimate 2: new velocity
    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

    // estimate 2: new residual force (and mass matrix if needed)
    ...

    // estimate 2: new acceleration
    M.solve(fr, a);

    // new velocity
    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

    // new residual force (and mass matrix if needed)
    ...

    // new acceleration
    M.solve(fr, a);
  }

Full example
^^^^^^^^^^^^

.. todo::

  Compile, run, and view instructions

:download:`main.cpp <./examples/Dynamics/ElasticLaminate-VelocityVerlet/main.cpp>`
:download:`CMakeLists.txt <./examples/Dynamics/ElasticLaminate-VelocityVerlet/CMakeLists.txt>`
:download:`plot.py <./examples/Dynamics/ElasticLaminate-VelocityVerlet/plot.py>`

Overdamped dynamics
===================

.. todo::

  Description and example

Forward Euler
-------------

.. code-block:: cpp

  // position and velocity
  xt::xtensor<double,2> u = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v = xt::zeros<double>({nnode, ndim});

  // time increments
  for ( ... )
  {
    // new displacement
    xt::noalias(u) = u + dt * v;

    // new velocity based on residual force
    ...
  }
