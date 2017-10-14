
.. _fem_examples_dynamic_diagonal-mass:

********************
Diagonal mass matrix
********************

'Theory'
========

[:download:`source: Tri3/complete.cpp <Tri3/complete.cpp>`]

[:download:`source: Tri3/lumped.cpp <Tri3/lumped.cpp>`]

[:download:`source: Tri3/diagonal.cpp <Tri3/diagonal.cpp>`]

[:download:`source: Quad4/complete.cpp <Quad4/complete.cpp>`]

[:download:`source: Quad4/lumped.cpp <Quad4/lumped.cpp>`]

[:download:`source: Quad4/diagonal.cpp <Quad4/diagonal.cpp>`]

.. note::

  Should you still want to use the full (sparse) matrix in conjunction with a solver, please consider using the Intel compile with the Intel MKL library (which includes the Pardiso solver). Consider `this example <https://github.com/tdegeus/cpp_examples/tree/master/eigen_mkl>`_.

Because dynamic computations need so many increments, it may be worth to cut down the computational cost of a simulation. The best way to start cutting is on the solver, which is often the most costly of the entire Finite Element Program. Remember that

.. math::

  \underline{\underline{M}}(\vec{x})\;
  \underline{\vec{a}}
  =
  \underline{\vec{t}}(\vec{x})
  -
  \underline{\vec{f}}(\vec{x})
  -
  \underline{\underline{H}}(\vec{x})\;
  \underline{\vec{v}}

If we make the mass matrix diagonal, :math:`\underline{\vec{a}}` can be obtained without solving a system. Instead we merely needed the component wise inverse of the diagonal terms (which are the only remaining non-zero terms).

There are essentially four possible definitions of the mass matrix:

1.  The *consistent* (or *complete*) mass matrix. Here the quadrature of the mass matrix is identical to that of the internal force in the system. It is the result of numerical integration of

    .. math::

      \underline{\underline{M}}(\vec{x})
      =
      \int\limits_\Omega
        \rho(\vec{x})\; \underline{N}(\vec{X})\; \underline{N}^\mathsf{T}(\vec{X}) \;
      \mathrm{d}\Omega

2.  The *lumped* mass matrix. In which

    .. math::

      M_{ii}^\mathrm{lumped} = \sum\limits_{j} M_{ij}

3.  The *lumped* mass matrix with *diagonal scaling*:

    .. math::

      M_{ii}^\mathrm{lumped} = c \; \sum\limits_{j} M_{ij}

    where the constant :math:`c` is set such that:

    .. math::

      \sum\limits_{i} M_{ii}^\mathrm{lumped}
      =
      \int\limits_\Omega
        \rho(\vec{x})\;
      \mathrm{d}\Omega

    A trick to accomplish this is to set

    .. math::

      c
      =
      \frac{
        \sum\limits_{i}\sum\limits_{j} M_{ij}
      }
      {
         \sum\limits_{i} M_{ii}
      }

    where the denominator is the trace of :math:`\underline{\underline{M}}`.

4.  Evaluating of :math:`\underline{\underline{M}}` using a quadrature rule involving only the nodal points and thus automatically yielding a diagonal matrix. It relies on the fact that when the quadrature-points coincide with the nodes, a local support is obtained (since the shape functions are zero in all the other nodes). The integration point volume is that part of the element's volume that can be associated to that node.

    For 2-D linear triangles, which have three nodes and normally one Gauss point in the middle of the element, the number of integration points is increased to three (whose positions correspond to that of the nodes), and their weight is set to 1/3 of the normal weight. Interestingly, for this case all possible definitions of the diagonalized mass matrix are the same.

    For 2-D quads only the position of the integration points changes to that of the nodes.

.. note:: References

  * `This answer on Computational Science StackExchange <https://scicomp.stackexchange.com/questions/19704/how-to-formulate-lumped-mass-matrix-in-fem>`_.




