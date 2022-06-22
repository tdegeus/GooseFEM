*************
Element::Hex8
*************

| :download:`GooseFEM/ElementHex8.h <../../include/GooseFEM/ElementHex8.h>`
| :download:`GooseFEM/ElementHex8.hpp <../../include/GooseFEM/ElementHex8.hpp>`

Element::Hex8::Quadrature
=========================

Element definition to numerically interpolate and integrate.

.. note::

  This function evaluates the shape function gradients upon construction, they are not recomputed upon evaluation. To evaluate them with respect to updated coordinates (e.g. to do updated Lagrange), use ".update_x(...)" to update the nodal coordinates and re-evaluate the shape function gradients and integration volumes.

.. note::

  By default integration is done using Gauss points. To use a different scheme one has to supply the position (in isoparametric coordinates) and weight of the integration points (their number is inferred from the input).

.. note::

  Most functions take the output as the last input-argument, as to write directly to a pre-allocated array, avoiding their re-allocation. All these functions have a wrapper that does the allocation for you (and thus returns the output rather than taking it as input). All function of this kind are indicated here with a *

Element::Hex8::Quadrature::update_x(...)
----------------------------------------

Update the nodal coordinates (elemvec: [nelem, nne, ndim]).

Element::Hex8::Quadrature::nelem()
----------------------------------

Number of elements.

Element::Hex8::Quadrature::nne()
--------------------------------

Number of nodes per element.

Element::Hex8::Quadrature::ndim()
---------------------------------

Number of dimensions.

Element::Hex8::Quadrature::nip()
--------------------------------

Number of integration points.

Element::Hex8::Quadrature::GradN()
----------------------------------

(Current) Shape function gradient (w.r.t. real coordinates): [nelem, nip, nne, ndim]

Element::Hex8::Quadrature::asTensor<...>(...)*
----------------------------------------------

Convert a 'qscalar' (scalar values stored per integration point) to a 'qtensor',
a tensor per integration point, with all tensor-components having the same value.
The template parameters allows you to specify the rank of the tensor.
From Python use the function that allocates data, and specify the rank as first
argument.

Element::Hex8::Quadrature::dV(...)
----------------------------------

(Current) Volume of each integration point (qscalar: [nelem, nip]).

Element::Hex8::Quadrature::gradN_vector(...)*
---------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \vec{\nabla} N_m \vec{u}_m

or in index notation

.. math::

  \varepsilon_{ij} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Hex8::Quadrature::gradN_vector_T(...)*
-----------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \left( \vec{\nabla} N_m \vec{u}_m \right)^T

or in index notation

.. math::

  \varepsilon_{ji} = \frac{\partial N_m}{\partial x_i} u_{mj}

Element::Hex8::Quadrature::symGradN_vector(...)*
------------------------------------------------

Implementation of

.. math::

  \bm{\varepsilon} = \tfrac{1}{2} \left(
    \vec{\nabla} N_m \vec{u}_m + \left( \vec{\nabla} N_m \vec{u}_m \right)^T
  \right)

Element::Hex8::Quadrature::int_N_scalar_NT_dV(...)*
---------------------------------------------------

Implementation of

.. math::

  M_{mn}
  =
  \int\limits_{\Omega^h} N_m \; \rho \; N_n \; \mathrm{d}\Omega^h
  \equiv
  \sum\limits_q \; N_m \; \rho \; N_n \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim]. This implies that all dimensions are the same.

Element::Hex8::Quadrature::int_gradN_dot_tensor2_dV(...)*
---------------------------------------------------------

Implementation of:

.. math::

  \vec{f}_m = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \bm{\sigma} \; \mathrm{d}\Omega^h

or in index notation

.. math::

  f_{mj} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} \sigma_{ij} \; \delta\Omega_q

Element::Hex8::Quadrature::int_gradN_dot_tensor4_dot_gradNT_dV(...)*
--------------------------------------------------------------------

Implementation of:

.. math::

  \bm{K}_{mn} = \int\limits_{\Omega^h} ( \vec{\nabla} N_m ) \cdot \mathbb{C} \cdot \vec{\nabla} N_n \; \mathrm{d}\Omega^h

or in index notation

.. math::

  K_{m+id, n+kd} = \sum\limits_q \; \frac{\partial N_m}{\partial x_i} C_{ijkl} \frac{\partial N_n}{\partial x_l} \; \delta\Omega_q

Note that the output is an "elemmat", which has shape [nelem, nne*ndim, nne*ndim].

Element::Hex8::Quadrature::AllocateQtensor<...>(...)
----------------------------------------------------

Allocate (and initialize) a 'qtensor' of a certain rank (template parameter).
From Python specify the rank as fist argument.

Element::Hex8::Quadrature::AllocateQscalar(...)
-----------------------------------------------

Shortcut for ``AllocateQtensor<0>(...)``.

Element::Hex8::Gauss
====================

Integration points according to exact integration using Gauss points.

Element::Hex8::Gauss::nip()
---------------------------

Returns the number of integration points.

Element::Hex8::Gauss::xi()
--------------------------

Returns the position of the integration points in isoparametric coordinates [nip, ndim] (with ndim = 3).

Element::Hex8::Gauss::w()
-------------------------

Returns the weights of the integration points [nip].

Element::Hex8::Nodal
====================

Integration points that coincide with the nodes (equally weight). This scheme can for example be used to obtain a diagonal mass matrix.

Element::Hex8::Nodal::nip()
---------------------------

Returns the number of integration points.

Element::Hex8::Nodal::xi()
--------------------------

Returns the position of the integration points in isoparametric coordinates [nip, ndim] (with ndim = 3).

Element::Hex8::Nodal::w()
-------------------------

Returns the weights of the integration points [nip].
