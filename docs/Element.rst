
*****************
GooseFEM::Element
*****************

Introduction
============

Provides routines to perform numerical quadrature. The philosophy is that the in- and output are multi-dimensional arrays (all of type ``GooseFEM::ArrD``) that contain:

+-----------+--------------------------------------------------------------+----------------------------------+
|  Alias    | Description                                                  | Shape                            |
+===========+==============================================================+==================================+
| "elemmat" | matrices stored per element                                  | [nelem, nne*ndim, nne*ndim]      |
+-----------+--------------------------------------------------------------+----------------------------------+
| "elemvec" | nodal vectors stored per element                             | [nelem, nne, ndim]               |
+-----------+--------------------------------------------------------------+----------------------------------+
| "qtensor" | tensors stored (as list) per integration point, per element  | [nelem, nip, #tensor-components] |
+-----------+--------------------------------------------------------------+----------------------------------+
| "qscalar" | scalars stored per integration point, per element            | [nelem, nip]                     |
+-----------+--------------------------------------------------------------+----------------------------------+

As a result of this choice, the routines here need to know nothing of your choice how to organise your data, and thus remain flexible on code-unspecific. Also, evaluation in parallel (with OpenMP) has been trivially implemented, without there being any pitfalls.

The different elements each have a class ``Quadrature`` that is used to take gradients and to integrate. The idea is that you supply the nodal coordinates as "elemvec" to the constructor, and if you wish customize the integration points and their weights. The shape functions and their gradients are immediately evaluated and stored for reuse.

A small example to get you started is to compute the displacement gradients at each integration point:

.. code-block:: cpp

  #include <GooseFEM/GooseFEM.h>
  #include <cppmat/cppmat.h>

  int main()
  {
    // some mesh
    GooseFEM::Mesh::Quad4::Regular mesh(3,3);

    // get relevant fields
    GooseFEM::MatD coor = mesh.coor();
    GooseFEM::MatS conn = mesh.conn();
    GooseFEM::MatS dofs = mesh.dofs();

    // define quadrature
    GooseFEM::Element::Quad4::Quadrature quad(GooseFEM::Element::asElementVector(conn,coor));

    // define some displacement field: simple shear
    // - zero-initialize displacements
    GooseFEM::MatD disp = GooseFEM::MatD::Zero(mesh.nnode(), mesh.ndim());
    // - shear stain
    double gamma = 0.1;
    // - update displacement field
    for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
      disp(n,0) = gamma * ( coor(n,1) - coor(0,1) );

    // compute the displacement gradient of each integration point of each element
    GooseFEM::ArrD Gradu = quad.gradN_vector_T(GooseFEM::Element::asElementVector(conn,disp));

    // view result
    // - tensor to interpret the tensor components of "Gradu": [nelem, nip, #tensor-components]
    cppmat::view::cartesian::tensor2<double,2> gradu;
    // - view for all integration points of all elements
    for ( size_t e = 0 ; e < mesh.nelem() ; ++e ) {
      for ( size_t k = 0 ; k < quad.nip() ; ++k ) {
        // -- interpret sub-matrix
        gradu.setMap(&Gradu(e,k));
        // -- print element and integration point number
        std::cout << "e = " << e << ", k = " << k << std::endl;
        // -- print the displacement gradient
        std::cout << std::setw(5) << std::setprecision(3) << gradu << std::endl;
      }
    }

    return 0;
  }

.. tip::

  :ref:`element_asElementVector` converts nodal coordinates and displacement to the corresponding "elemvec". In addition:

  *   :ref:`vector`: switch between "nodevec", "elemvec", "dofval", ...

.. tip::

  To take the gradients and integral with respect to updated coordinates (i.e. to do updated Lagrange), use the ``.update_x(...)`` method to update the nodal coordinates and re-evaluate the shape function gradients and integration volumes.

.. note::

  All routines that take or return a "qtensor" have been templated such that you can supply a type to interpret it. For example:





.. note::

  The code and headers for the different elements are quite similar. They have been kept as parallel implementations to allow flexible adaption. One can inspect or deploy changes easily using an editor that highlight the differences between files.


GooseFEM::Element::Quad4
========================

No description yet, please consult the code.

GooseFEM::Element::Quad4::Quadrature
------------------------------------

No description yet, please consult the code.

GooseFEM::Element::Quad4::Quadrature::gradN_vector
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No description yet, please consult the code.

GooseFEM::Element::Quad4::Nodal
-------------------------------

No description yet, please consult the code.

GooseFEM::Element::Quad4::Gauss
-------------------------------

No description yet, please consult the code.

General routines
================

.. _element_asElementVector:

GooseFEM::Element::asElementVector
----------------------------------

No description yet, please consult the code.

