
********
Dynamics
********

Time discretisation
===================

Verlet
------

.. code-block:: cpp

  while ( ... )
  {
    // variables & history

    xt::xtensor<double,2> u;
    xt::xtensor<double,1> V;
    xt::xtensor<double,1> A;
    xt::xtensor<double,1> V_n = g.dofs_v();
    xt::xtensor<double,1> A_n = g.dofs_a();

    // new displacement

    xt::noalias(u) = g.u() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a();

    g.set_u(u);

    // new acceleration

    xt::noalias(A) = g.solve_A();

    g.set_a(A);

    // new velocity

    xt::noalias(V) = V_n + .5 * dt * ( A_n + A );

    g.set_v(V);
  }

Velocity Verlet
---------------

.. code-block:: cpp

  while ( ... )
  {
    // variables & history

    xt::xtensor<double,2> u;
    xt::xtensor<double,1> V;
    xt::xtensor<double,1> A;
    xt::xtensor<double,1> V_n = g.dofs_v();
    xt::xtensor<double,1> A_n = g.dofs_a();

    // new displacement

    xt::noalias(u) = g.u() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a();

    g.set_u(u);

    // estimate new velocity

    xt::noalias(V) = V_n + dt * A_n;

    g.set_v(V);

    xt::noalias(A) = g.solve_A();

    xt::noalias(V) = V_n + .5 * dt * ( A_n + A );

    g.set_v(V);

    // new velocity

    xt::noalias(A) = g.solve_A();

    xt::noalias(V) = V_n + .5 * dt * ( A_n + A );

    g.set_v(V);

    // new acceleration

    xt::noalias(A) = g.solve_A();

    g.set_a(A);
  }

Forward Euler
-------------

.. code-block:: cpp

  while ( ... )
  {
    xt::noalias(U) = g.dofs_u() + dt * g.solve_V();

    g.set_u(U);
  }
