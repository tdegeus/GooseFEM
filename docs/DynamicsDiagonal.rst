
****************************
GooseFEM::Dynamics::Diagonal
****************************

[:download:`source: DynamicsDiagonal.h <../src/GooseFEM/DynamicsDiagonal.h>`, :download:`source: DynamicsDiagonal.cpp <../src/GooseFEM/DynamicsDiagonal.cpp>`]

Overview
========

Principle
---------

The philosophy is to provide some structure to efficiently run a finite element simulation which remains customizable. Even more customization can be obtained by copy/pasting the source and modifying it to your need. The idea that is followed involves a hierarchy of three classes, whereby a class that is higher in hierarchy writes to some field of the class that is lower in hierarchy, runs a function, and reads from some field. In general:

*   **Discretized system** (``GooseFEM::Dynamics::Diagonal::Periodic``, ``GooseFEM::Dynamics::Diagonal::SemiPeriodic``).

    *   Defines the discretized system.
    *   Writes element positions, displacements, and velocity of all elements to the *element definition*.
    *   Assembles the diagonal (inverse) mass matrix, the displacement dependent forces, the velocity dependent forces, and the diagonal damping matrix from the element arrays computed in *element definition*.
    *   Provides time integrators.

*   **Element definition** (``GooseFEM::Dynamics::Diagonal::SmallStrain::Qaud4``)

    Provides the element arrays by performing numerical quadrature. At the integration point the strain and strain-rate are computed and constitutive response is probed from the quadrature point definition.

*   **Quadrature point definition**

    This class is not provided, and should be provided by the user.

Example
-------

A simple example is:

[:download:`source: examples/DynamicsDiagonalPeriodic/laminate/no_damping/Verlet/main.cpp <examples/DynamicsDiagonalPeriodic/laminate/no_damping/Verlet/main.cpp>`]

.. code-block:: cpp

  #include <Eigen/Eigen>
  #include <cppmat/cppmat.h>
  #include <GooseFEM/GooseFEM.h>
  #include <GooseMaterial/AmorphousSolid/SmallStrain/Elastic/Cartesian2d.h>

  // -------------------------------------------------------------------

  using     MatS = GooseFEM::MatS;
  using     MatD = GooseFEM::MatD;
  using     ColD = GooseFEM::ColD;

  using     T2   = cppmat::cartesian2d::tensor2 <double>;
  using     T2s  = cppmat::cartesian2d::tensor2s<double>;

  namespace GM   = GooseMaterial::AmorphousSolid::SmallStrain::Elastic::Cartesian2d;

  // ===================================================================

  class Quadrature
  {
  public:
    T2s eps, epsdot, sig;

    size_t nhard;
    GM::Material hard, soft;

    Quadrature(size_t nhard);

    double density             (size_t e, size_t k);
    void   stressStrain        (size_t e, size_t k);
    void   stressStrainRate    (size_t e, size_t k);
    void   stressStrainPost    (size_t e, size_t k);
    void   stressStrainRatePost(size_t e, size_t k);
  };

  // -------------------------------------------------------------------

  Quadrature::Quadrature(size_t _nhard)
  {
    nhard    = _nhard;
    hard     = GM::Material(100.,10.);
    soft     = GM::Material(100., 1.);
  }

  // -------------------------------------------------------------------

  double Quadrature::density(size_t elem, size_t k, double V)
  {
    return 1.0;
  }
  // -------------------------------------------------------------------

  void Quadrature::stressStrain(size_t elem, size_t k, double V)
  {
    if ( elem < nhard ) sig = hard.stress(eps);
    else                sig = soft.stress(eps);
  }
  // -------------------------------------------------------------------

  void Quadrature::stressStrainRate(size_t elem, size_t k, double V)
  {
  }
  // -------------------------------------------------------------------

  void Quadrature::stressStrainPost(size_t elem, size_t k, double V)
  {
    Vbar += V;

    if ( elem < nhard ) Ebar += hard.energy(eps) * V;
    else                Ebar += soft.energy(eps) * V;
  }

  // -------------------------------------------------------------------

  void Quadrature::stressStrainRatePost(size_t elem, size_t k, double V)
  {
  }

  // ===================================================================

  int main()
  {
    // class which provides the mesh
    GooseFEM::Mesh::Quad4::Regular mesh(40,40,1.);

    // class which provides the constitutive response at each quadrature point
    auto  quadrature = std::make_shared<Quadrature>(40*40/4);

    // class which provides the response of each element
    using Elem = GooseFEM::Dynamics::Diagonal::SmallStrain::Quad4<Quadrature>;
    auto  elem = std::make_shared<Elem>(quadrature);

    // class which provides the system and an increment
    GooseFEM::Dynamics::Diagonal::Periodic<Elem> sim(
      elem,
      mesh.coor(),
      mesh.conn(),
      mesh.dofsPeriodic(),
      1.e-2,
      0.0
    );

    // loop over increments
    for ( ... )
    {
      // - set displacement of fixed DOFs
      ...

      // - compute time increment
      sim.Verlet();

      // - post-process
      quadrature->Ebar = 0.0;
      quadrature->Vbar = 0.0;

      sim.post();

      ...
    }

    return 0;
  }

Pseudo-code
-----------

What is happening inside ``Verlet`` is evaluating the forces (and the mass matrix), and updating the displacements by solving the system. In pseudo-code:

*   Mass matrix:

    .. code-block:: python

      sim.computeMinv():
      {
        for e in elements:

          sim->elem->xe(i,j) = ...
          sim->elem->ue(i,j) = ...

          sim->elem->computeM(e):
          {
            for k in integration-points:

              sim->elem->M(...,...) += ... * sim->elem->quad->density(e,k,V)
          }

          M(...) += sim->elem->M(i,i)
      }

*   Displacement dependent force:

    .. code-block:: python

      sim.computeFu():
      {
        for e in elements:

          sim->elem->xe(i,j) = ...
          sim->elem->ue(i,j) = ...

          sim->elem->computeFu(e):
          {
            for k in integration-points:

              sim->elem->quad->eps(i,j) = ...

              sim->elem->quad->stressStrain(e,k,V)

              sim->elem->fu(...) += ... * sim->elem->quad->sig(i,j)
          }

          Fu(...) += sim->elem->fu(i)
      }

*   Velocity dependent force:

    .. code-block:: python

      sim.computeFv():
      {
        for e in elements:

          sim->elem->xe(i,j) = ...
          sim->elem->ue(i,j) = ...
          sim->elem->ve(i,j) = ...

          sim->elem->computeFv(e):
          {
            for k in integration-points:

              sim->elem->quad->epsdot(i,j) = ...

              sim->elem->quad->stressStrainRate(e,k,V)

              sim->elem->fv(...) += ... * sim->elem->quad->sig(i,j)
          }

          Fv(...) += sim->elem->fu(i)
      }

Signature
---------

From this it is clear that:

*   ``GooseFEM::Dynamics::Diagonal::Periodic`` requires the following minimal signature from ``GooseFEM::Dynamics::Diagonal::SmallStrain::Qaud4``:

    .. code-block:: cpp

      class Element
      {
      public:
        matrix M;                    // should have operator(i,j)
        column fu, fv;               // should have operator(i)
        matrix xe, ue, ve;           // should have operator(i,j)

        void computeM (size_t elem); // mass matrix                     <- quad->density
        void computeFu(size_t elem); // displacement dependent forces   <- quad->stressStrain
        void computeFv(size_t elem); // displacement dependent forces   <- quad->stressStrainRate
        void post     (size_t elem); // post-process                    <- quad->stressStrain(Rate)
      }

*   ``GooseFEM::Dynamics::Diagonal::SmallStrain::Qaud4`` requires the minimal signature from ``Quadrature``

    .. code-block:: cpp

      class Quadrature
      {
      public:
        tensor eps, epsdot, sig;     // should have operator(i,j)

        double density             (size_t elem, size_t k, double V);
        void   stressStrain        (size_t elem, size_t k, double V);
        void   stressStrainRate    (size_t elem, size_t k, double V);
        void   stressStrainPost    (size_t elem, size_t k, double V);
        void   stressStrainRatePost(size_t elem, size_t k, double V);
      }

