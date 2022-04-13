**************
Element::Quad4
**************

| :download:`GooseFEM/ElementQuad4.h <../../include/GooseFEM/ElementQuad4.h>`
| :download:`GooseFEM/ElementQuad4.hpp <../../include/GooseFEM/ElementQuad4.hpp>`
| :download:`GooseFEM/ElementQuad4Planar.h <../../include/GooseFEM/ElementQuad4Planar.h>`
| :download:`GooseFEM/ElementQuad4Planar.hpp <../../include/GooseFEM/ElementQuad4Planar.hpp>`
| :download:`GooseFEM/ElementQuad4Axisymmetric.h <../../include/GooseFEM/ElementQuad4Axisymmetric.h>`
| :download:`GooseFEM/ElementQuad4Axisymmetric.hpp <../../include/GooseFEM/ElementQuad4Axisymmetric.hpp>`

Element::Quad4::Quadrature
==========================

Element definition to numerically interpolate and integrate.

.. note::

  This function evaluates the shape function gradients upon construction, they are not recomputed upon evaluation. To evaluate them with respect to updated coordinates (e.g. to do updated Lagrange), use ".update_x(...)" to update the nodal coordinates and re-evaluate the shape function gradients and integration volumes.

.. note::

  By default integration is done using Gauss points. To use a different scheme one has to supply the position (in isoparametric coordinates) and weight of the integration points (their number is inferred from the input).

.. note::

  Most functions take the output as the last input-argument, as to write directly to a pre-allocated array, avoiding their re-allocation. All these functions have a wrapper that does the allocation for you (and thus returns the output rather than taking it as input). All function of this kind are indicated here with a *

Element::Quad4::Quadrature::update_x(...)
-----------------------------------------

Update the nodal coordinates (elemvec: [nelem, nne, ndim]).

Element::Quad4::Quadrature::nelem()
-----------------------------------

Number of elements.

Element::Quad4::Quadrature::nne()
---------------------------------

Number of nodes per element.

Element::Quad4::Quadrature::ndim()
----------------------------------

Number of dimensions.

Element::Quad4::Quadrature::nip()
---------------------------------

Number of integration points.

Element::Quad4::Quadrature::GradN()
-----------------------------------

(Current) Shape function gradient (w.r.t. real coordinates): [nelem, nip, nne, ndim]

Element::Quad4::Quadrature::asTensor<...>(...)*
-----------------------------------------------

Convert a 'qscalar' (scalar values stored per integration point) to a 'qtensor',
a tensor per integration point, with all tensor-components having the same value.
The template parameters allows you to specify the rank of the tensor.
From Python use the function that allocates data, and specify the rank as first
argument.

Element::Quad4::Quadrature::dV(...)
-----------------------------------

(Current) Volume of each integration point (qscalar: [nelem, nip]).

Element::Quad4::Quadrature::gradN_vector(...)*
----------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \vec{\nabla} N_m \vec{u}_m

or in index notation

.. math::

  \varepsilon_{ij} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Quad4::Quadrature::gradN_vector_T(...)*
------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \left( \vec{\nabla} N_m \vec{u}_m \right)^T

or in index notation

.. math::

  \varepsilon_{ji} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Quad4::Quadrature::symGradN_vector(...)*
-------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \tfrac{1}{2} \left(
    \vec{\nabla} N_m \vec{u}_m + \left( \vec{\nabla} N_m \vec{u}_m \right)^T
  \right)

Element::Quad4::Quadrature::int_N_scalar_NT_dV(...)*
----------------------------------------------------

Implementation of

.. math::

  M_{mn}
  =
  \int\limits_{\Omega^h} N_m \; \rho \; N_n \; \mathrm{d}\Omega^h
  \equiv
  \sum\limits_q \; N_m \; \rho \; N_n \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim]. This implies that all dimensions are the same.

Element::Quad4::Quadrature::int_gradN_dot_tensor2_dV(...)*
----------------------------------------------------------

Implementation of:

.. math::

  \vec{f}_m = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \bm{\sigma} \; \mathrm{d}\Omega^h

or in index notation

.. math::

  f_{mj} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} \sigma_{ij} \; \delta\Omega_q

Element::Quad4::Quadrature::int_gradN_dot_tensor4_dot_gradNT_dV(...)*
---------------------------------------------------------------------

Implementation of:

.. math::

  \bm{K}_{mn} = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \mathbb{C} \cdot \vec{\nabla} N_n \; \mathrm{d}\Omega^h

or in index notation

.. math::

  K_{m+id, n+kd} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} C_{ijkl} \frac{\partial N_n}{\partial x_l} \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim].

Element::Quad4::QuadraturePlanar
================================

Element definition to numerically interpolate and integrate under a planar assumption. This implies that all the tensors are 3-d, but that the third dimension is ignored by all functions (although for output these components are zero-initialised).

.. note::

  This function evaluates the shape function gradients upon construction, they are not recomputed upon evaluation. To evaluate them with respect to updated coordinates (e.g. to do updated Lagrange), use ".update_x(...)" to update the nodal coordinates and re-evaluate the shape function gradients and integration volumes.

.. note::

  By default integration is done using Gauss points. To use a different scheme one has to supply the position (in isoparametric coordinates) and weight of the integration points (their number is inferred from the input).

.. note::

  Most functions take the output as the last input-argument, as to write directly to a pre-allocated array, avoiding their re-allocation. All these functions have a wrapper that does the allocation for you (and thus returns the output rather than taking it as input). All function of this kind are indicated here with a *

Element::Quad4::QuadraturePlanar::update_x(...)
-----------------------------------------------

Update the nodal coordinates (elemvec: [nelem, nne, ndim]).

Element::Quad4::QuadraturePlanar::nelem()
-----------------------------------------

Number of elements.

Element::Quad4::QuadraturePlanar::nne()
---------------------------------------

Number of nodes per element.

Element::Quad4::QuadraturePlanar::ndim()
----------------------------------------

Number of dimensions.

Element::Quad4::QuadraturePlanar::nip()
---------------------------------------

Number of integration points.

Element::Quad4::QuadraturePlanar::GradN()
-----------------------------------------

(Current) Shape function gradient (w.r.t. real coordinates): [nelem, nip, nne, ndim]

Element::Quad4::QuadraturePlanar::dV(...)*
------------------------------------------

(Current) Volume of each integration point (qscalar: [nelem, nip]). An overload is available to get the same result as a tensor per integration point (qtensor: [nelem, nip, tdim, tdim]) with all tensor-components having the same value.

Element::Quad4::QuadraturePlanar::gradN_vector(...)*
----------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \vec{\nabla} N_m \vec{u}_m

or in index notation

.. math::

  \varepsilon_{ij} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Quad4::QuadraturePlanar::gradN_vector_T(...)*
------------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \left( \vec{\nabla} N_m \vec{u}_m \right)^T

or in index notation

.. math::

  \varepsilon_{ji} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Quad4::QuadraturePlanar::symGradN_vector(...)*
-------------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \tfrac{1}{2} \left(
    \vec{\nabla} N_m \vec{u}_m + \left( \vec{\nabla} N_m \vec{u}_m \right)^T
  \right)

Element::Quad4::QuadraturePlanar::int_N_scalar_NT_dV(...)*
----------------------------------------------------------

Implementation of

.. math::

  M_{mn}
  =
  \int\limits_{\Omega^h} N_m \; \rho \; N_n \; \mathrm{d}\Omega^h
  \equiv
  \sum\limits_q \; N_m \; \rho \; N_n \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim]. This implies that all dimensions are the same.

Element::Quad4::QuadraturePlanar::int_gradN_dot_tensor2_dV(...)*
----------------------------------------------------------------

Implementation of:

.. math::

  \vec{f}_m = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \bm{\sigma} \; \mathrm{d}\Omega^h

or in index notation

.. math::

  f_{mj} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} \sigma_{ij} \; \delta\Omega_q

Element::Quad4::QuadraturePlanar::int_gradN_dot_tensor4_dot_gradNT_dV(...)*
---------------------------------------------------------------------------

Implementation of:

.. math::

  \bm{K}_{mn} = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \mathbb{C} \cdot \vec{\nabla} N_n \; \mathrm{d}\Omega^h

or in index notation

.. math::

  K_{m+id, n+kd} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} C_{ijkl} \frac{\partial N_n}{\partial x_l} \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim].

Element::Quad4::QuadratureAxisymmetric
======================================

Element definition to numerically interpolate and integrate in an axisymmetric cylindrical coordinate system. This implies that all tensors (stress, strain, stiffness, ...) are fully three dimensional, but the discretisation is two-dimensional.

.. note::

  This function evaluates the shape function gradients upon construction, they are not recomputed upon evaluation. To evaluate them with respect to updated coordinates (e.g. to do updated Lagrange), use ".update_x(...)" to update the nodal coordinates and re-evaluate the shape function gradients and integration volumes.

.. note::

  By default integration is done using Gauss points. To use a different scheme one has to supply the position (in isoparametric coordinates) and weight of the integration points (their number is inferred from the input).

.. note::

  Most functions take the output as the last input-argument, as to write directly to a pre-allocated array, avoiding their re-allocation. All these functions have a wrapper that does the allocation for you (and thus returns the output rather than taking it as input). All function of this kind are indicated here with a *

Element::Quad4::QuadratureAxisymmetric::update_x(...)
-----------------------------------------------------

Update the nodal coordinates (elemvec: [nelem, nne, ndim]).

Element::Quad4::QuadratureAxisymmetric::nelem()
-----------------------------------------------

Number of elements.

Element::Quad4::QuadratureAxisymmetric::nne()
---------------------------------------------

Number of nodes per element.

Element::Quad4::QuadratureAxisymmetric::ndim()
----------------------------------------------

Number of dimensions.

Element::Quad4::QuadratureAxisymmetric::nip()
---------------------------------------------

Number of integration points.

Element::Quad4::QuadratureAxisymmetric::GradN()
-----------------------------------------------

(Current) Shape function gradient (w.r.t. real coordinates): [nelem, nip, nne, ndim]

Element::Quad4::QuadratureAxisymmetric::dV(...)*
------------------------------------------------

(Current) Volume of each integration point (qscalar: [nelem, nip]). An overload is available to get the same result as a tensor per integration point (qtensor: [nelem, nip, tdim, tdim]) with all tensor-components having the same value.

Element::Quad4::QuadratureAxisymmetric::gradN_vector(...)*
----------------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \vec{\nabla} N_m \vec{u}_m

or in index notation

.. math::

  \varepsilon_{ij} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Quad4::QuadratureAxisymmetric::gradN_vector_T(...)*
------------------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \left( \vec{\nabla} N_m \vec{u}_m \right)^T

or in index notation

.. math::

  \varepsilon_{ji} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Quad4::QuadratureAxisymmetric::symGradN_vector(...)*
-------------------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \tfrac{1}{2} \left(
    \vec{\nabla} N_m \vec{u}_m + \left( \vec{\nabla} N_m \vec{u}_m \right)^T
  \right)

Element::Quad4::QuadratureAxisymmetric::int_N_scalar_NT_dV(...)*
----------------------------------------------------------------

Implementation of

.. math::

  M_{mn}
  =
  \int\limits_{\Omega^h} N_m \; \rho \; N_n \; \mathrm{d}\Omega^h
  \equiv
  \sum\limits_q \; N_m \; \rho \; N_n \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim]. This implies that all dimensions are the same.

Element::Quad4::QuadratureAxisymmetric::int_gradN_dot_tensor2_dV(...)*
----------------------------------------------------------------------

Implementation of:

.. math::

  \vec{f}_m = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \bm{\sigma} \; \mathrm{d}\Omega^h

or in index notation

.. math::

  f_{mj} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} \sigma_{ij} \; \delta\Omega_q

Element::Quad4::QuadratureAxisymmetric::int_gradN_dot_tensor4_dot_gradNT_dV(...)*
---------------------------------------------------------------------------------

Implementation of:

.. math::

  \bm{K}_{mn} = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \mathbb{C} \cdot \vec{\nabla} N_n \; \mathrm{d}\Omega^h

or in index notation

.. math::

  K_{m+id, n+kd} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} C_{ijkl} \frac{\partial N_n}{\partial x_l} \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim].

Element::Quad4::Quadrature::AllocateQtensor<...>(...)
-----------------------------------------------------

Allocate (and initialize) a 'qtensor' of a certain rank (template parameter).
From Python specify the rank as fist argument.

Element::Quad4::Quadrature::AllocateQscalar(...)
------------------------------------------------

Shortcut for ``AllocateQtensor<0>(...)``.

Element::Quad4::Gauss
=====================

Integration points according to exact integration using Gauss points.

Element::Quad4::Gauss::nip()
----------------------------

Returns the number of integration points.

Element::Quad4::Gauss::xi()
---------------------------

Returns the position of the integration points in isoparametric coordinates [nip, ndim] (with ndim = 3).

Element::Quad4::Gauss::w()
--------------------------

Returns the weights of the integration points [nip].

Element::Quad4::Nodal
=====================

Integration points that coincide with the nodes (equally weight). This scheme can for example be used to obtain a diagonal mass matrix.

Element::Quad4::Nodal::nip()
----------------------------

Returns the number of integration points.

Element::Quad4::Nodal::xi()
---------------------------

Returns the position of the integration points in isoparametric coordinates [nip, ndim] (with ndim = 3).

Element::Quad4::Nodal::w()
--------------------------

Returns the weights of the integration points [nip].

Element::Quad4::MidPoint
========================

Single integration point in the middle of the element.

Element::Quad4::MidPoint::nip()
-------------------------------

Returns the number of integration points.

Element::Quad4::MidPoint::xi()
------------------------------

Returns the position of the integration points in isoparametric coordinates [nip, ndim] (with ndim = 3).

Element::Quad4::MidPoint::w()
-----------------------------

Returns the weights of the integration points [nip].
