
**************
Linear elastic
**************

Consider a uniform linear elastic bar that is extended by a
uniform fixed displacement on both sides.
This problem can be modelled and discretised using symmetry as show below.
In this example we will furthermore assume that the bar is sufficiently thick
in the out-of-plane direction to be modelled using two-dimensional plane strain.

.. image:: statics/FixedDisplacements_LinearElastic/problem-sketch.svg
    :width: 300px
    :align: center

|

Below, an example is described line-by-line.
For simplicity the material model is taken from
`GMatElastic <https://www.github.com/tdegeus/GMatElastic>`_, as this example is about FEM,
not about constitutive modelling.
The full example can be downloaded run as follows:

.. tabs::

    .. tab:: C++

        :download:`CMakeLists.txt <statics/FixedDisplacements_LinearElastic/CMakeLists.txt>`
        :download:`example.cpp <statics/FixedDisplacements_LinearElastic/example.cpp>`
        :download:`plot.py <statics/FixedDisplacements_LinearElastic/plot.py>`

        To compile we can use for example *CMake*:

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/CMakeLists.txt
            :language: cmake
            :lines: 1-13

        Whereby we notice the use of *GMatElastic* as described,
        but also *HighFive* illustrate a way of storing data and then plotting it with Python.

        Proceed for example as follows:

        1.  Get the prerequisites::

                conda install -c conda-forge goosefem gmatelastic highfive python goosempl

        2.  Compile::

                cmake -Bbuild
                cd build
                cmake --build .

        3.  Run::

                ./example

        4.  Plot::

                python ../plot.py

    .. tab:: Python

        | :download:`example.py <statics/FixedDisplacements_LinearElastic/example.py>`

        In this example we illustrate some basic plotting.

        Proceed for example as follows:

        1.  Get the prerequisites::

                conda install -c conda-forge python-goosefem python-gmatelastic goosempl

        2.  Run (and plot)::

                python example.py

Include library
===============

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 1-4
            :emphasize-lines: 2-3

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 1-9
            :emphasize-lines: 4

The first step is to include the header-only library.
Some dependencies are included for convenience.

Define mesh
===========

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 11-30
            :emphasize-lines: 2

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 14-33
            :emphasize-lines: 2

A mesh is defined using *GooseFEM*.
As observed the *mesh* is a class that has methods to extract the relevant information
such as the nodal coordinates (*coor*), the connectivity (*conn*),
the degrees-of-freedom per node (*dofs*) and several node-sets that
will be used to impose the sketched boundary conditions
(*nodesLft*, *nodesRgt*, *nodesTop*, *nodesBot*).

Note that:

*   The connectivity (*conn*) contains information of which nodes, in which order,
    belong to which element.

*   The degrees-of-freedom per node (*dofs*) contains information of how a nodal vector
    (a vector stored per node) can be transformed to a list of degrees-of-freedom as used
    in the linear system (although this can be mostly done automatically as we will see below).

.. seealso::

  * :ref:`conventions_terminology`
  * Details: :ref:`MeshQuad4`

Define partitioning
===================

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 35-40

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 38-43

We will reorder such that degrees-of-freedom are ordered such that

.. math::

    \texttt{u} =
    \begin{bmatrix}
        \texttt{u}_u \\
        \texttt{u}_p
    \end{bmatrix}

where the subscript :math:`u` and :math:`p` respectively denote
*Unknown* and *Prescribed*
degrees-of-freedom.
To achieve this we start by collecting all prescribed degrees-of-freedom in *iip*.

(Avoid) Book-keeping
====================

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 45
            :emphasize-lines: 1

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 49
            :emphasize-lines: 1

To switch between the three of *GooseFEM*'s data-representations,
an instance of the *Vector* class is used.
This instance, *vector*, will enable us to switch between a vector field (e.g. the displacement)

1.  collected per node,
2.  collected per degree-of-freedom, or
3.  collected per element.

.. note::

    The *Vector* class collects most, if not all, the burden of book-keeping.
    It is thus here that *conn*, *dofs*, and *iip* are used. In particular,

    *   'nodevec' :math:`\leftrightarrow` 'dofval' using *dofs* and *iip*,
    *   'nodevec' :math:`\leftrightarrow` 'elemvec' using *conn*.

    By contrast, most of *GooseFEM*'s other methods receive the relevant representation,
    and consequently require no problem specific knowledge.
    They thus do not have to supplied with *conn*, *dofs*, or *iip*.

.. seealso::

    *   :ref:`conventions_vector`
    *   :ref:`conventions_storage`
    *   Details: :ref:`Vector`

System matrix
=============

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 48-49
            :emphasize-lines: 1

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 52-53
            :emphasize-lines: 1

We now also allocate the system/stiffness system (stored as sparse matrix).
Like vector, it can accept and return different vector representations,
in addition to the ability to assemble from element system matrices.

In addition we allocate the accompanying sparse solver,
that we will use to solve a linear system of equations.
Note that the solver-class takes care of factorising only when needed
(when the matrix has been changed).

.. tabs::

    .. tab:: C++

        .. note::

            Here, the default solver is used (which is the default template, hence the "<>").
            To use other solvers see: :ref:`linear_solver`.

.. seealso::

    *   :ref:`conventions_matrix`
    *   Details: :ref:`Matrix`

Allocate nodal & element vectors
================================

To avoid repeated memory allocation,
it is advised to pre-allocate some data array and reuse them.
We allocate:

*   *fint*: nodal internal forces
*   *fres*: nodal residual forces

and the following arrays (tensors per element), that eliminate the connectivity from the equation,
and allow a generic API:

*   *ue*: displacement per element
*   *fe*: force per element (strictly speaking *ue* could be reused)
*   *Ke*: tangent matrix per element.

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 51-58

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 55-62

.. tip::

    To allocate nodal vectors or tensors per element, use the convenience functions:

    .. tabs::

        .. tab:: C++

            .. code-block:: cpp

                // nodal vectors ("fint", "fext", "fext", "disp", or "coor")
                auto shape = vector.shape_nodevec(); // get shape
                auto variable = vector.allocate_nodevec(); // allocate
                auto variable = vector.allocate_elemvec(0.0); // allocate & (zero-)initialise

                // vector per element ("ue" or "fe")
                auto shape = vector.shape_elemvec(); // shape
                auto variable = vector.allocate_elemvec(); // allocate
                auto variable = vector.allocate_elemvec(0.0); // allocate & (zero-)initialise

                // matrix per element ("Ke")
                auto shape = vector.shape_elemmat(); // shape
                auto variable = vector.allocate_elemmat(); // allocate
                auto variable = vector.allocate_elemmat(0.0); // allocate & (zero-)initialise

        .. tab:: Python

            .. code-block:: python

                # nodal vectors ("fint", "fext", "fext", "disp", or "coor")
                variable = np.zeros(vector.shape_nodevec())

                # vector per element ("ue" or "fe")
                variable = np.zeros(vector.shape_elemvec())

                # matrix per element ("Ke")
                variable = np.zeros(vector.shape_elemmat())

.. warning::

    Upsizing (e.g. *disp* :math:`\rightarrow` *ue*) can be done uniquely,
    but downsizing (e.g. *fe* :math:`\rightarrow` *fint*) can be done in more than one way,
    see :ref:`conventions_vector_conversion`.
    We will get back to this point below.


Element definition
==================

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 64-65

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 68-69

At this moment the interpolation and quadrature is allocated.
The shape functions and integration points (that can be customised) are stored in this class.
As observed, no further information is needed than the number of elements and
the nodal coordinates per element.
Both are contained in the output of ``vector.AsElement(coor)``, which is an 'elemvec' of
shape "[nelem, nne, ndim]".
This illustrates that problem specific book-keeping is isolated to the main program,
using *Vector* as tool.

.. note::

    The shape-functions are computed when constructing this class,
    they are not recomputed when evaluating them.
    One can recompute them if the nodal coordinates change using ".update_x(...)", however,
    this is only relevant in a large deformation setting.

.. seealso::

    *   :ref:`conventions_vector`
    *   :ref:`conventions_storage`
    *   Details: :ref:`Vector`
    *   Details: :ref:`ElementQuad4`

Material definition
===================

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 68

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 72

We now define a uniform linear elastic material,
using an external library that is tuned to *GooseFEM*.
This material library will translate a strain tensor per integration point to a stress tensor
per integration point and a stiffness tensor per integration point.

.. seealso::

    Material libraries tuned to *GooseFEM* include:

    *   `GMatElastic <https:://www.github.com/tdegeus/GMatElastic>`__
    *   `GMatElastoPlastic <https:://www.github.com/tdegeus/GMatElastoPlastic>`__
    *   `GMatElastoPlasticFiniteStrainSimo <https:://www.github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo>`__
    *   `GMatElastoPlasticQPot <https:://www.github.com/tdegeus/GMatElastoPlasticQPot>`__
    *   `GMatElastoPlasticQPot3d <https:://www.github.com/tdegeus/GMatElastoPlasticQPot3d>`__
    *   `GMatNonLinearElastic <https:://www.github.com/tdegeus/GMatNonLinearElastic>`__

    But other libraries can also be easily used with (simple) wrappers.

Allocate integration point tensors
==================================

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 71-74

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 75-77

We will need a few tensors per integration.
Like before, we can choose to allocate them to avoid repeated memory allocation.
In particular, we allocate the strain, stress, and stiffness tensors per integration point.
Note that these tensors are 3-d while our problem was 2-d.
This is thanks to the plane strain assumption,
and the element definition that ignores all third-dimension components.

.. note::

    To allocate integration point tensors, use the convenience functions:

    .. tabs::

        .. tab:: C++

            .. code-block:: cpp

                auto shape = quad.shape_qtensor<rank>(); // shape
                auto variable = quad.allocate_qtensor<rank>(); // allocate
                auto variable = quad.allocate_qtensor<rank>(0.0); // allocate & (zero-)initialise

        .. tab:: Python

            .. code-block:: python

                variable = np.zeros(quad.shape_qtensor(rank))

Compute strain
==============

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 79-81
            :emphasize-lines: 2

    .. tab:: C++ (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.cpp
            :language: cpp
            :lines: 65

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 83-84
            :emphasize-lines: 2

    .. tab:: Python (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.py
            :language: python
            :lines: 69

The strain per integration point is now computed using the current nodal displacements
(stored as 'elemvec' in *ue*) and the gradient of the shape functions.

.. note::

    *ue* is the output of ``vector.asElement(disp, ue)``.
    Using this syntax re-allocation of *ue* is avoided.
    If this optimisation is irrelevant for you problem,
    please use the same function, but starting with a **capital letter**.

Compute stress and tangent
==========================

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 83-85

    .. tab:: C++ (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.cpp
            :language: cpp
            :lines: 68-70

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 87-89

    .. tab:: Python (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.py
            :language: python
            :lines: 72-74

The stress and stiffness tensors are now computed for each integration point
(completely independently) using the external material model.

Assemble system
===============

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 87-94

    .. tab:: C++ (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.cpp
            :language: cpp
            :lines: 72-76

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 91-97

    .. tab:: Python (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.py
            :language: python
            :lines: 76-80

The stress stored per integration point (*Sig*) is now converted to
nodal internal force vectors stored per element (*fe*).
Using *vector* this 'elemvec' representation is then converted of a
'nodevec' representation in *fint*.
Likewise, the stiffness tensor stored for integration point (*C*) are converted
to system matrices stored per element ('elemmat') and finally assembled to
the global stiffness matrix.

.. warning::

    Please note that downsizing
    (*fe* :math:`\rightarrow` *fint* and *Ke* :math:`\rightarrow` *K*) can be done in two ways,
    and that "assemble..." is the right function here as it adds entries that occur
    more than once.
    In contrast "as..." would not result in what we want here.

Solve
=====

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 95-106

    .. tab:: C++ (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.cpp
            :language: cpp
            :lines: 78-88

    .. tab:: C++ (manual)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/manual_partition.cpp
            :language: cpp
            :lines: 78-95

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 99-109

    .. tab:: Python (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.py
            :language: python
            :lines: 82-92

    .. tab:: Python (manual)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/manual_partition.py
            :language: python
            :lines: 82-102

We now prescribe the displacement of the Prescribed degrees-of-freedom directly
in the nodal displacements *disp* and compute the residual force.
This is follows by partitioning and solving, all done internally in the *MatrixPartitioned* class.
As an example, the same operation with manual book-keeping is included.

Post-process
============

Strain and stress
-----------------

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 111-114

    .. tab:: C++ (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.cpp
            :language: cpp
            :lines: 94-96

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 115-118

    .. tab:: Python (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.py
            :language: python
            :lines: 98-100

The strain and stress per integration point are recomputed for post-processing.

Residual force
--------------

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 115-127

    .. tab:: C++ (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.cpp
            :language: cpp
            :lines: 98-108

    .. tab:: C++ (manual)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/manual_partition.cpp
            :language: cpp
            :lines: 105-118

    .. tab:: Python

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 120-131

    .. tab:: Python (realloc)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example_realloc.py
            :language: python
            :lines: 102-112

    .. tab:: Python (manual_partition)

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/manual_partition.py
            :language: python
            :lines: 112-125

We convince ourselves that the solution is indeed in mechanical equilibrium.

Plot
----

.. tabs::

    .. tab:: C++

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.cpp
            :language: cpp
            :lines: 128-138

        Finally we store some fields for plotting using
        :download:`plot.py <statics/FixedDisplacements_LinearElastic/plot.py>`.

    .. tab:: Python

        Let's extract the average stress per element:

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 133-135

        And plot it on a deformed mesh:

        .. literalinclude:: statics/FixedDisplacements_LinearElastic/example.py
            :language: python
            :lines: 147-186
