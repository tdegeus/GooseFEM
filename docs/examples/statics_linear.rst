
**************
Linear elastic
**************

Consider a uniform linear elastic bar that is extended by a uniform fixed displacement on both sides. This problem can be modelled and discretised using symmetry as show below. In this example we will furthermore assume that the bar is sufficiently thick in the out-of-plane direction to be modelled using two-dimensional plane strain.

.. image:: statics/FixedDisplacements_LinearElastic/problem-sketch.svg
  :width: 300px
  :align: center

|

Below an example is described line-by-line. The full example can be downloaded:

| :download:`main.cpp <statics/FixedDisplacements_LinearElastic/example/main.cpp>`
| :download:`CMakeLists.txt <statics/FixedDisplacements_LinearElastic/example/CMakeLists.txt>`
| :download:`plot.py <statics/FixedDisplacements_LinearElastic/example/plot.py>`

.. todo::

  Compile and run instructions.

.. note::

  This example is also available using the Python interface (:download:`main.py <statics/FixedDisplacements_LinearElastic/example/main.py>`). Compared to the C++ API, the Python API requires more data-allocation, in particular for the functions "AsElement" and "AssembleNode". See: :ref:`conventions_allocation`.

Include library
===============

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 1-4
   :emphasize-lines: 1-2

The first step is to include the header-only library. Note that for this example we also make use of a material model (`GMatElastic <https://www.github.com/tdegeus/GMatElastic>`_) and a library to write (and read) HDF5 files (`HighFive <https://www.github.com/BlueBrain/HighFive>`_).

Define mesh
===========

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 11-28
   :emphasize-lines: 2

A mesh is defined using GooseFEM. As observed the "mesh" is a class that has methods to extract the relevant information such as the nodal coordinates ("coor"), the connectivity ("conn"), the degrees-of-freedom per node ("dofs") and several node-sets that will be used to impose the sketched boundary conditions ("nodesLeft", "nodesRight", "nodesTop", "nodesBottom").

Note that:

* The connectivity ("conn") contains information of which nodes, in which order, belong to which element.
* The degrees-of-freedom per node ("dofs") contains information of how a nodal vector (a vector stored per node) can be transformed to a list of degrees-of-freedom as used in the linear system (although this can be mostly done automatically as we will see below).

.. seealso::

  * :ref:`conventions_terminology`
  * Details: :ref:`MeshQuad4`

Define partitioning
===================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 33-38

We will reorder such that degrees-of-freedom are ordered such that

.. math::

  \texttt{u} =
  \begin{bmatrix}
    \texttt{u}_u \\
    \texttt{u}_p
  \end{bmatrix}

where the subscript :math:`u` and :math:`p` respectively denote *Unknown* and *Prescribed* degrees-of-freedom. To achieve this we start by collecting all prescribed degrees-of-freedom in "iip".

(Avoid) Book-keeping
====================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 44
   :emphasize-lines: 1

To switch between the three of GooseFEM's data-representations, an instance of the "Vector" class is used. This instance, "vector", will enable us to switch between a vector field (e.g. the displacement)

1. collected per node,
2. collected per degree-of-freedom, or
3. collected per element.

.. note::

  The "Vector" class collects most, if not all, the burden of book-keeping. It is thus here that "conn", "dofs", and "iip" are used. In particular,

  * 'nodevec' :math:`\leftrightarrow` 'dofval' using "dofs" and "iip",
  * 'nodevec' :math:`\leftrightarrow` 'elemvec' using "conn".

  By contrast, most of GooseFEM's other methods receive the relevant representation, and consequently require no problem specific knowledge. They thus do not have to supplied with "conn", "dofs", or "iip".

.. seealso::

  * :ref:`conventions_vector`
  * :ref:`conventions_storage`
  * Details: :ref:`Vector`

System matrix
=============

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 47
   :emphasize-lines: 1

We now also allocate the system/stiffness system (stored as sparse matrix). Like vector, it can accept and return different vector representations, in addition to the ability to assemble from element system matrices.

.. note::

  Here, the default solver is used (which is the default template, hence the "<>"). To use other solvers see: :ref:`linear_solver`.

.. seealso::

  * :ref:`conventions_matrix`
  * Details: :ref:`Matrix`

Allocate nodal vectors
======================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 50-53

* "disp": nodal displacements
* "fint": nodal internal forces
* "fext": nodal external forces
* "fres": nodal residual forces

Allocate element vectors
========================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 56-58

* "ue": displacement
* "fe": force
* "Ke": tangent matrix

.. warning::

  Upsizing (e.g. "disp" :math:`\rightarrow` "ue") can be done uniquely, but downsizing (e.g. "fe" :math:`\rightarrow` "fint") can be done in more than one way, see :ref:`conventions_vector_conversion`. We will get back to this point below.

Element definition
==================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 64-65

At this moment the interpolation and quadrature is allocated. The shape functions and integration points (that can be customised) are stored in this class. As observed, no further information is needed than the number of elements and the nodal coordinates per element. Both are contained in the output of "vector.AsElement(coor)", which is an 'elemvec' of shape "[nelem, nne, ndim]". This illustrates that problem specific book-keeping is isolated to the main program, using "Vector" as tool.

.. note::

  The shape functions are computed when constructing this class, they are not recomputed when evaluating them. One can recompute them if the nodal coordinates change using ".update_x(...)", however, this is only relevant in a large deformation setting.

.. seealso::

  * :ref:`conventions_vector`
  * :ref:`conventions_storage`
  * Details: :ref:`Vector`
  * Details: :ref:`ElementQuad4`

Material definition
===================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 68

We now define a uniform linear elastic material, using an external library that is tuned to GooseFEM. This material library will translate a strain tensor per integration point to a stress tensor per integration point and a stiffness tensor per integration point.

.. seealso::

  Material libraries tuned to GooseFEM include:

  * `GMatElastic <https:://www.github.com/tdegeus/GMatElastic>`__
  * `GMatElastoPlastic <https:://www.github.com/tdegeus/GMatElastoPlastic>`__
  * `GMatElastoPlasticFiniteStrainSimo <https:://www.github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo>`__
  * `GMatElastoPlasticQPot <https:://www.github.com/tdegeus/GMatElastoPlasticQPot>`__
  * `GMatElastoPlasticQPot3d <https:://www.github.com/tdegeus/GMatElastoPlasticQPot3d>`__
  * `GMatNonLinearElastic <https:://www.github.com/tdegeus/GMatNonLinearElastic>`__

  But other libraries can also be easily used with (simple) wrappers.

Integration point tensors
=========================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 71-74

These strain, stress, and stiffness tensors per integration point are now allocated. Note that these tensors are 3-d while our problem was 2-d. This is thanks to the plane strain assumption, and the element definition that ignores all third-dimension components.

Compute strain
==============

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 80-81

The strain per integration point is now computed using the current nodal displacements (stored as 'elemvec' in "ue") and the gradient of the shape functions.

.. note::

  "ue" is the output of "vector.asElement(disp, ue)". Using this syntax re-allocation of "ue" is avoided. If this optimisation is irrelevant for you problem (or if you are using the Python interface), please use the same function, but starting with a capital:

  .. code-block:: cpp

    ue = vector.AsElement(disp);

  Note that this allows the one-liner

  .. code-block:: cpp

    Eps = elem.SymGradN_vector(vector.AElement(disp));

Compute stress and tangent
==========================

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 84

The stress and stiffness tensors are now computed for each integration point (completely independently) using the external material model.

.. note::

  "Sig" and "C" are the output variables that were preallocated in the main.

Assemble system
===============

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 86-92

The stress stored per integration point ("Sig") is now converted to nodal internal force vectors stored per element ("fe"). Using "vector" this 'elemvec' representation is then converted of a 'nodevec' representation in "fint". Likewise, the stiffness tensor stored for integration point ("C") are converted to system matrices stored per element ('elemmat') and finally assembled to the global stiffness matrix.

.. warning::

  Please note that downsizing ("fe" :math:`\rightarrow` "fint" and "Ke" :math:`\rightarrow` "K") can be done in two ways, and that "assemble..." is the right function here as it adds entries that occur more than once. In contrast "as..." would not result in what we want here.

.. note::

  Once more, "fe", "fint", and "Ke" are output variables. Less efficient, but shorter, is:

  .. code-block:: cpp

    // internal force
    fint = vector.AssembleNode(elem.Int_gradN_dot_tensor2_dV(Sig));

    // stiffness matrix
    K.assemble(elem.Int_gradN_dot_tensor4_dot_gradNT_dV(C));

Solve
=====

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 94-104

We now prescribe the displacement of the Prescribed degrees-of-freedom directly in the nodal displacements "disp" and compute the residual force. This is follows by partitioning and solving, all done internally in the "MatrixPartitioned" class.

Post-process
============

Strain and stress
-----------------

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 110-112

The strain and stress per integration point are recomputed for post-processing.

Residual force
--------------

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 114-125

We convince ourselves that the solution is indeed in mechanical equilibrium.

Store & plot
------------

.. literalinclude:: statics/FixedDisplacements_LinearElastic/example/main.cpp
   :language: cpp
   :lines: 127-136

Finally we store some fields for plotting using :download:`plot.py <statics/FixedDisplacements_LinearElastic/example/plot.py>`.

Manual partitioning
===================

To verify how partitioning and solving is done internally using the "MatrixPartitioned" class, the same example is provided where partitioning is done manually:

| :download:`main.cpp <statics/FixedDisplacements_LinearElastic/manual_partition/main.cpp>`
| :download:`CMakeLists.txt <statics/FixedDisplacements_LinearElastic/manual_partition/CMakeLists.txt>`
| :download:`plot.py <statics/FixedDisplacements_LinearElastic/manual_partition/plot.py>`
| :download:`main.py <statics/FixedDisplacements_LinearElastic/manual_partition/main.py>`
