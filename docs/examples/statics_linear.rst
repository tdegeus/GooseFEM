
**************
Linear elastic
**************

.. note::

    For readability the examples are only presented on the Python API.
    Converting them to C++ is trivial: the API is (almost) identical,
    and also the API of *xtensor* is almost identical to that of *NumPy*.

Consider a uniform linear elastic bar that is extended by a
uniform fixed displacement on both sides.
This problem can be modelled and discretised using symmetry as show below.
In this example we will furthermore assume that the bar is sufficiently thick
in the out-of-plane direction to be modelled using two-dimensional plane strain.

| :download:`sketch.svg <statics/fixed-displacement_elastic_sketch.svg>`

Below, an example is described line-by-line.
For simplicity the material model is taken from
`GMatElastic <https://www.github.com/tdegeus/GMatElastic>`_, as this example is about FEM,
not about constitutive modelling.
The full example can be downloaded run as follows:

| :download:`example.py <statics/fixed-displacement_elastic.py>`

In this example we illustrate some basic plotting.

Proceed for example as follows:

1.  Get the prerequisites::

        conda install -c conda-forge python-goosefem python-gmatelastic goosempl

2.  Run (and plot)::

        python fixed-displacement_elastic.py --plot

Include library
===============

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 1-6
    :emphasize-lines: 5

Define mesh
===========

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 11-24
    :emphasize-lines: 2

A mesh is defined using *GooseFEM*.
As observed the *mesh* is a class that has methods to extract the relevant information
such as the nodal coordinates (*coor*), the connectivity (*conn*),
the degrees-of-freedom per node (*dofs*) and several node-sets that
will be used to impose the sketched boundary conditions
(*nodesRightEdge*, *nodesTopEdge*, *nodesLeftEdge*, *nodesBottomEdge*).

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

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 27-34

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

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 40
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

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 43-33
    :emphasize-lines: 1

We now also allocate the system/stiffness system (stored as sparse matrix).
Like vector, it can accept and return different vector representations,
in addition to the ability to assemble from element system matrices.

In addition we allocate the accompanying sparse solver,
that we will use to solve a linear system of equations.
Note that the solver-class takes care of factorising only when needed
(when the matrix has been changed).

.. seealso::

    *   :ref:`conventions_matrix`
    *   Details: :ref:`Matrix`

Element definition
==================

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 47-48

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

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 53

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

Compute strain
==============

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 59-61

The strain per integration point is now computed using the current nodal displacements
(stored as 'elemvec' in *ue*) and the gradient of the shape functions.

.. note::

    The displacement per element *ue* eliminate the connectivity from the equation
    for a very significant part of the API.
    This makes things much easier to control and much more modular.

.. warning::

    Upsizing (e.g. *disp* :math:`\rightarrow` *ue*) can be done uniquely,
    but downsizing (e.g. *fe* :math:`\rightarrow` *fint*) can be done in more than one way,
    see :ref:`conventions_vector_conversion`.
    We will get back to this point below.

.. note::

    A function like

    .. literalinclude:: statics/fixed-displacement_elastic.py
        :language: python
        :lines: 59

    starting with a **capital letter**, allocates the output data.
    However, if that data was already allocated, it can be modified in-place with the
    corresponding function starting with a **lowercase letter**, e.g.

    .. literalinclude:: statics/fixed-displacement_elastic.py
        :language: python
        :lines: 87

.. tip::

    To allocate nodal vectors or tensors per element, use the convenience functions:

    .. code-block:: python

        # nodal vectors ("fint", "fext", "fext", "disp", or "coor")
        variable = np.zeros(vector.shape_nodevec())

        # vector per element ("ue" or "fe")
        variable = np.zeros(vector.shape_elemvec())

        # matrix per element ("Ke")
        variable = np.zeros(vector.shape_elemmat())

Assemble system
===============

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 64-69

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

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 71-81

We now compute the residual forces (trivially zero in this case, but kept for pedagogical reasons).
We then prescribe the displacement of the Prescribed degrees-of-freedom directly
in the nodal displacements *disp* and compute the residual force.
This is follows by partitioning and solving, all done internally in the *MatrixPartitioned* class.
As an example, the same operation with manual book-keeping is included in

| :download:`example.py <statics/fixed-displacement_elastic_manual-partition.py>`

Post-process
============

Strain and stress
-----------------

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 87-89

The strain and stress per integration point are recomputed for post-processing.

Residual force
--------------

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 91-102

We convince ourselves that the solution is indeed in mechanical equilibrium.

Plot
----

Let's extract the average (equivalent) stress per element:

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 120-122

And plot an equivalent stress on a deformed mesh:

.. literalinclude:: statics/fixed-displacement_elastic.py
    :language: python
    :lines: 126-127
