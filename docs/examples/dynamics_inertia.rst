
*******
Inertia
*******

Time discretisation: Verlet
===========================

.. code-block:: cpp

    // position, velocity, acceleration (and history: last increment)
    xt::xtensor<double,2> u   = xt::zeros<double>({nnode, ndim});
    xt::xtensor<double,2> v   = xt::zeros<double>({nnode, ndim});
    xt::xtensor<double,2> a   = xt::zeros<double>({nnode, ndim});
    xt::xtensor<double,2> v_n = xt::zeros<double>({nnode, ndim});
    xt::xtensor<double,2> a_n = xt::zeros<double>({nnode, ndim});

    // residual force
    xt::xtensor<double,2> fres = xt::zeros<double>({nnode, ndim});

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
        xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt, 2.0) * a;

        // new residual force (and mass matrix if needed)
        ...

        // new acceleration
        M.solve(fres, a);

        // new velocity
        xt::noalias(v) = v_n + 0.5 * dt * (a_n + a);
    }

Example
=======

.. todo::

  Compile, run, and view instructions:

  .. code-block:: bash

    mkdir build
    cd build
    make
    ./main
    python3 plot.py

| :download:`main.cpp <dynamics/Elastic-Verlet/main.cpp>`
| :download:`CMakeLists.txt <dynamics/Elastic-Verlet/CMakeLists.txt>`
| :download:`plot.py <dynamics/Elastic-Verlet/plot.py>`

