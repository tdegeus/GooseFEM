**********
Overdamped
**********

Forward Euler
=============

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
