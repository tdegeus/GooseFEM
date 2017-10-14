
.. _fem_examples_small-strain_nonlinear_dense:

**********************************
Non-linear statics -- small strain
**********************************

.. contents:: **Contents**
  :local:
  :depth: 2
  :backlinks: top

Here we extend the example from :ref:`fem_examples_small-strain_linear_dense` to a non-linear constitutive response, which is however still subjected to a **small deformations** assumption. We treat the constitutive response as a black-box here:

.. math::

  \bm{\varepsilon}
  \;
  \rightarrow
  \;
  \bm{\sigma}, \mathbb{K}

where it must be emphasized that :math:`\mathbb{K}` symmetrizes. To understand more about the constitutive response, please consult the `documentation of GooseSolid <https://github.com/tdegeus/GooseSolid/blob/master/docs/LinearElastic/NonLinearElastic.pdf>`_

Mixed boundary conditions
=========================

[:download:`source: fixedbnd.cpp <fixedbnd.cpp>`]

In summary, our iterative update reads

.. math::

  \underline{\underline{\bm{K}}}_{(i)} \cdot \delta \underline{\vec{x}}
  =
  \underline{\vec{t}}
  -
  \underline{\vec{f}}_{(i)}

with

.. math::

  \underline{\underline{\bm{K}}}_{(i)}
  =
  \int\limits_{\Omega^h_0}
    \big[\, \vec{\nabla}_0 \underline{N} \,\big]
    \cdot
    \bm{K}\big(\vec{x}_{(i)}\big)
    \cdot
    \big[\, \vec{\nabla}_0 \underline{N} \,\big]^\mathsf{T} \;
  \mathrm{d}\Omega

and

.. math::

  \underline{\vec{f}}_{(i)}
  =
  \int\limits_{\Omega^h_0}
    \big[\, \vec{\nabla}_0 \underline{N} \,\big]
    \cdot
    \bm{\sigma}\big(\vec{x}_{(i)}\big) \;
  \mathrm{d}\Omega

We will use this to iteratively update

.. math::

  \underline{\vec{x}}_{(i+1)} = \underline{\vec{x}}_{(i)} + \delta \underline{\vec{x}}

From an 'initial guess'

.. math::

  \underline{\vec{x}}_{(0)} = \underline{\vec{0}}

.. note::

  This is a bit of particular case. We need this iteration to get things going, however typically

  .. math::

    \underline{\vec{t}}
    -
    \underline{\vec{f}}_{(0)}
    =
    \underline{\vec{0}}

  One should not be fooled, this does not mean to an equilibrium has been obtained.

Like before, we introduce DOFs to make the system scalar. We then need to partition the system to deal with the prescribed displacement-components:

.. math::

  \begin{bmatrix}
  \underline{\underline{K}}_{uu}^{(i)} && \underline{\underline{K}}_{up}^{(i)} \\
  \underline{\underline{K}}_{pu}^{(i)} && \underline{\underline{K}}_{pp}^{(i)}
  \end{bmatrix}
  \cdot
  \begin{bmatrix}
  \delta \underline{x}_{u}^{(i)} \\
  \delta \underline{x}_{p}^{(i)}
  \end{bmatrix}
  =
  \begin{bmatrix}
  \underline{t}_{u}^{(i)} \\
  \underline{t}_{p}^{(i)}
  \end{bmatrix}
  -
  \begin{bmatrix}
  \underline{f}_{u}^{(i)} \\
  \underline{f}_{p}^{(i)}
  \end{bmatrix}

From which the update of the unknown DOFs as follows

.. math::

  \delta \underline{x}_{u}^{(i)}
  =
  \left[
    \underline{\underline{K}}_{uu}^{(i)}
  \right]^{-1}
  \left[
    \underline{t}_{u}^{(i)} -
    \underline{f}_{u}^{(i)} -
    \underline{\underline{K}}_{up}^{(i)} \delta \underline{x}_{p}^{(i)}
  \right]

There a very important concept hidden here. Because we prescribe :math:`\delta \underline{x}_{p}` directly to the correct value -- which we do not iteratively update -- we need to set

.. math::

  \delta \underline{x}_{p}^{(i)}
  =
  \begin{cases}
    \delta \underline{x}_{p} \quad & \mathrm{if}\, i = 0 \\
    0                        \quad & \mathrm{otherwise}
  \end{cases}

(Which means that one does not have to compute the product :math:`\underline{\underline{K}}_{up}^{(i)} \delta \underline{x}_{p}^{(i)}` for :math:`i > 0`, as it will be zero.)

.. note:: **Reaction forces**

  To obtain the reaction forces on the prescribed DOFs simply use that

  .. math::

    \underline{t}_{p}^{(i)} = \underline{f}_{p}^{(i)}

  Which should be evaluated once convergence has been reached (before that, this has no meaning).

Periodic problem
================

[:download:`source: periodic.cpp <periodic.cpp>`]

In some ways the periodic example is even easier than the example above. As discussed previously :math:`\delta \underline{\vec{x}}` only leads to periodic fluctuations. In summary we begin by setting

.. math::

  \underline{\vec{x}}_{(0)}
  =
  \big[\, \bar{\bm{F}} - \bm{I} \,\big]
  \cdot
  \big[\, \underline{\vec{X}} - \underline{\vec{X}}_\mathrm{ref} \,\big]

the update

.. math::

  \underline{\vec{x}}_{(i+1)} = \underline{\vec{x}}_{(i)} + \delta \underline{\vec{x}}

then only affects the fluctuations, and not the average. Also we do not need to worry about the above discussion on the prescribed displacements since they are always zero and traction free -- as they merely suppress rigid body modes.

.. note::

  The periodic implementation with macroscopic DOFs is identically extended to the non-linear case as above. Here one does have care about properly adding the prescribed displacements.
