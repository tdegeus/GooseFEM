

*   `GooseMaterial <https://github.com/tdegeus/GooseMaterial>`_

    Provides the constitutive response (and optionally the constitutive tangent) of several materials.


*   `cppmat <https://github.com/tdegeus/cppmat>`_

    Provides tensor classes and operations. (The number of tensor operations are limited in the main program, and even non-standard, but this library is crucial to compute the material response implemented in `GooseMaterial <https://github.com/tdegeus/GooseMaterial>`_.)


*   `Eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_

    A linear algebra library. As you will notice, Eigen plays an important role in GooseFEM, and glues everything together since in the end the Finite Element Method is just a way to cast a problem into a set linear or linearized equations. Most of the efficiency of the final program will depend on the efficiency of the implementation of the linear algebra. In several examples we will simplify the structure by using dense matrices together with a simple solver which solves the resulting linear system. In reality one should always use sparse matrices combined with a more efficient solver. As you will notice, many examples need only few changes to be transformed in a production code.

.. note:: **Compilation**

  Unless otherwise mentioned, the examples can be compiled as follows. Provided that ``pkg-config`` is set-up correctly one can use

  .. code-block:: bash

    clang++ `pkg-config --cflags Eigen3 cppmat GooseMaterial GooseFEM` -std=c++14 -o example example_name.cpp

  (whereby ``clang++`` can be replaced by for example ``g++``). If one does not want to use ``pkg-config``, one has to specify ``-I/path/to/library`` for each of the libraries.

  For further development it is strongly advised to include the options ``-Wpedantic -Wall`` to get on top of mistakes. Once the code is ready, one should compile with optimization (``-O3``) and without assertions (``-DNDEBUG``). The `Eigen3 documentation <http://eigen.tuxfamily.org/index.php?title=FAQ#How_can_I_enable_vectorization.3F>`_ further recommends the option ``-march=native`` to enable vectorization optimized for you architecture.




.. note::

  This is a good point to study some examples:

  *   :ref:`fem_examples_small-strain_linear_dense`

      We slowly work up to an iterative scheme starting from a linear problem, written, however, in such a way that the step towards a non-linear problem is small.

  *   :ref:`fem_examples_small-strain_nonlinear_dense`

      Here we employ Newton-Raphson to solve the non-linear equilibrium equation. It is easy to see that once the above examples have been understood this step is indeed trivial.

See: :ref:`fem_examples_dynamic_diagonal-mass`.






.. math::

  \rho\, \ddot{\vec{x}}
  =
  \vec{\nabla} \cdot
  \bm{\sigma}(\vec{x})
  +
  \eta\, \nabla^2\dot{\vec{x}}
  \qquad
  \vec{x} \in \Omega

where :math:`\rho` is the density and :math:`\eta` the viscosity (a.k.a. the damping coefficient). The first and second time derivative of the position :math:`\vec{x}` are respectively the velocity :math:`\vec{v} = \dot{\vec{x}}` and the acceleration :math:`\vec{a} = \ddot{\vec{x}}`.

We can generalize this as follows (which will also simplify our proceedings below)

.. math::

  \rho(\vec{x})\, \ddot{\vec{x}}
  =
  \vec{\nabla} \cdot
  \big[\, \bm{\sigma}(\vec{x}) + \bm{\sigma}_{\eta}(\vec{\dot{x}} ) \,\big]
  \qquad
  \vec{x} \in \Omega

.. note::

  To retrieve the original form

  .. math::

    \bm{\sigma}_{\eta} = \eta\; \vec{\nabla} \dot{\vec{x}}

  But, we can now also use other expressions. For example, the damping equivalent of linear elasticity:

  .. math::

    \bm{\sigma}_{\eta} (\vec{x}) = \mathbb{C}_{\eta} (\vec{x}) : \dot{\bm{\varepsilon}} (\vec{x})

  with

  .. math::

    \mathbb{C}_{\eta} (\vec{x})
    =
    \kappa (\vec{x}) \bm{I} \otimes \bm{I}
    +
    2 \gamma (\vec{x}) \mathbb{I}_d

  where :math:`\kappa` is the bulk viscosity while :math:`\gamma` is the shear viscosity. Furthermore

  .. math::

    \dot{\bm{\varepsilon}} (\vec{x})
    =
    \tfrac{1}{2} \big[\, \vec{\nabla} \dot{\vec{x}} + [\, \vec{\nabla} \dot{\vec{x}} \,]^T \,\big]

  Our original form is retrieved when :math:`\kappa = \tfrac{2}{3} \gamma`, both are independent of :math:`\vec{x}`, and :math:`\dot{\vec{x}}` possesses the necessary symmetries.




Like before, we will solve this equation in a weak sense

.. math::

  \int\limits_\Omega
    \rho(\vec{x})\; \vec{\phi}(\vec{X}) \cdot \ddot{\vec{x}} \;
  \mathrm{d}\Omega
  =
  \int\limits_\Omega
    \vec{\phi}(\vec{X})
    \cdot
    \Big[\,
      \vec{\nabla}
      \cdot
      \big[\, \bm{\sigma}(\vec{x}) + \bm{\sigma}_{\eta}(\vec{\dot{x}} ) \,\big]
    \,\Big] \;
  \mathrm{d}\Omega
  \qquad
  \forall \; \vec{\phi}(\vec{X}) \in \mathbb{R}^d

Integration by parts results in

.. math::

  \int\limits_\Omega
    \rho(\vec{x})\; \vec{\phi}(\vec{X}) \cdot \ddot{\vec{x}} \;
  \mathrm{d}\Omega
  =
  \int\limits_\Gamma
    \vec{\phi}(\vec{X}) \cdot \big[\, \vec{t}(\vec{x}) + \vec{t}_{\eta}(\vec{x}) \,\big] \;
  \mathrm{d}\Gamma
  -
  \int\limits_\Omega
    \big[\, \vec{\nabla} \vec{\phi}(\vec{X}) \,\big]
    :
    \big[\, \bm{\sigma}(\vec{x}) + \bm{\sigma}_{\eta}(\dot{\vec{x}}) \,\big] \;
  \mathrm{d}\Omega
  \qquad
  \forall \; \vec{\phi}(\vec{X}) \in \mathbb{R}^d

Which we will discretize as before:

.. math::

  \underline{\vec{\phi}}^\mathsf{T} \cdot
  \int\limits_\Omega
    \rho(\vec{x})\; \underline{N}(\vec{X})\; \underline{N}^\mathsf{T}(\vec{X}) \;
  \mathrm{d}\Omega \;
  \underline{\ddot{\vec{x}}}
  =
  \underline{\vec{\phi}}^\mathsf{T} \cdot
  \int\limits_\Gamma
    \underline{N}(\vec{X})\; \big[\, \vec{t}(\vec{x}) + \vec{t}_{\eta}(\vec{x}) \,\big] \;
  \mathrm{d}\Gamma
  -
  \underline{\vec{\phi}}^\mathsf{T} \cdot
  \int\limits_\Omega
    \big[\, \vec{\nabla} \underline{N}(\vec{X}) \,\big]
    :
    \big[\, \bm{\sigma}(\vec{x}) + \bm{\sigma}_{\eta}(\dot{\vec{x}}) \,\big] \;
  \mathrm{d}\Omega
  \qquad
  \forall \; \underline{\vec{\phi}} \in \mathbb{R}^d_n

Which is independent of the test functions, hence:

.. math::

  \int\limits_\Omega
    \rho(\vec{x})\; \underline{N}(\vec{X})\; \underline{N}^\mathsf{T}(\vec{X}) \;
  \mathrm{d}\Omega \;
  \underline{\ddot{\vec{x}}}
  =
  \int\limits_\Gamma
    \underline{N}(\vec{X})\; \big[\, \vec{t}(\vec{x}) + \vec{t}_{\eta}(\vec{x}) \,\big] \;
  \mathrm{d}\Gamma
  -
  \int\limits_\Omega
    \big[\, \vec{\nabla} \underline{N}(\vec{X}) \,\big]
    :
    \big[\, \bm{\sigma}(\vec{x}) + \bm{\sigma}_{\eta}(\dot{\vec{x}}) \,\big] \;
  \mathrm{d}\Omega

Which we can denote as follows

.. math::

  \underline{\underline{M}}(\vec{x})\; \underline{\ddot{\vec{x}}}
  =
  \underline{\vec{t}}(\vec{x})
  +
  \underline{\vec{t}}_{\eta}(\vec{x})
  -
  \underline{\vec{f}}(\vec{x})
  -
  \underline{\vec{f}}_{\eta}(\vec{x})

whereby we have introduced:

*   *Mass matrix*

    .. math::

      \underline{\underline{M}}(\vec{x})
      =
      \int\limits_\Omega
        \rho(\vec{x})\; \underline{N}(\vec{X})\; \underline{N}^\mathsf{T}(\vec{X}) \;
      \mathrm{d}\Omega

*   *Boundary tractions*

    .. math::

      \underline{\vec{t}}(\vec{x})
      =
      \int\limits_\Gamma
        \underline{N}(\vec{X})\; \vec{t}(\vec{x}) \;
      \mathrm{d}\Gamma
      \qquad
      \mathrm{and}
      \qquad
      \underline{\vec{t}}_{\eta}(\vec{x})
      =
      \int\limits_\Gamma
        \underline{N}(\vec{X})\; \vec{t}_{\eta}(\vec{x}) \;
      \mathrm{d}\Gamma

*   *Internal forces*

    .. math::

      \underline{\vec{f}}(\vec{x})
      =
      \int\limits_\Omega
        \big[\, \vec{\nabla} \underline{N}(\vec{X}) \,\big] : \bm{\sigma}(\vec{x}) \;
      \mathrm{d}\Omega
      \qquad
      \mathrm{and}
      \qquad
      \underline{\vec{f}}(\vec{x})
      =
      \int\limits_\Omega
        \big[\, \vec{\nabla} \underline{N}(\vec{X}) \,\big] : \bm{\sigma}_{\eta}(\dot{\vec{x}}) \;
      \mathrm{d}\Omega

.. note::

  In many problems it makes sense to assume the mass matrix constant, as any change of volume results in an equivalent change of the density, i.e.

  .. math::

    \int\limits_{\Omega}
      \rho(\vec{x})
    \;\mathrm{d}\Omega
    =
    \int\limits_{\Omega_0}
      \rho(\vec{X})
    \;\mathrm{d}\Omega_0

  This results in the following expression for the mass matrix:

  .. math::

    \underline{\underline{M}}(\vec{X})
    =
    \int\limits_{\Omega_0}
      \rho(\vec{X})\; \underline{N}(\vec{X})\; \underline{N}^\mathsf{T}(\vec{X}) \;
    \mathrm{d}\Omega_0
    =
    \mathrm{constant}

Time discretization
-------------------
.. note:: References

  `Syllabus of the course "Computational Physics (PY 502)" by Anders Sandvik, Department of Physics, Boston University <http://physics.bu.edu/py502/syllabus.pdf>`_.











.. note::

  The details depend on the element type. Several standard elements types are implemented in `GooseFEM <https://github.com/tdegeus/GooseFEM>`_.

