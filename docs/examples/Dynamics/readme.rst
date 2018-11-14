
********
Dynamics
********

Time discretisation
===================

Verlet
------

.. code-block:: cpp

  xt::xtensor<double,2> u    = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v    = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> a    = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v_n  = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> a_n  = xt::zeros<double>({nnode, ndim});

  xt::xtensor<double,2> fres = xt::zeros<double>({nnode, ndim});

  GooseFEM::MatrixDiagonal M(...);

  ...

  while ( ... )
  {
    // history

    xt::noalias(v_n) = v;
    xt::noalias(a_n) = a;

    // new displacement

    xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt,2.) * a;

    ...

    // new acceleration

    M.solve(fres, a);

    // new velocity

    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

  }

Velocity Verlet
---------------

.. code-block:: cpp

  while ( ... )
  {
    // variables & history

    xt::noalias(v_n) = v;
    xt::noalias(a_n) = a;

    // new displacement

    xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt,2.) * a;

    ...

    // estimate new velocity

    xt::noalias(v) = v_n + dt * a_n;

    ...

    M.solve(fres, a);

    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

    ...

    // new velocity

    M.solve(fres, a);

    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

    ...

    // new acceleration

    M.solve(fres, a);

  }

Forward Euler
-------------

.. code-block:: cpp

  xt::xtensor<double,2> u = xt::zeros<double>({nnode, ndim});
  xt::xtensor<double,2> v = xt::zeros<double>({nnode, ndim});

  ...

  while ( ... )
  {
    xt::noalias(u) = u + dt * v;

    ...
  }
