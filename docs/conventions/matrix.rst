.. _conventions_matrix:

Matrix representation
=====================

Element system matrix
---------------------

The element system matrix collects individual system matrices as a
multi-dimensional array of shape
:math:`\left[ n_\text{elements} \times n_\text{nodes-per-element} n_\text{dim} \times n_\text{nodes-per-element} n_\text{dim} \right]`.
An element system matrix

.. code-block:: python

    Ke = K[e,:,:]

obeys the following convention:

.. math::

    \underline{f} = \underline{\underline{K}} \underline{u}

where

.. math::

    f_{(n + i d)} \equiv f_i^{(n)}

with :math:`n` the node number, :math:`i` the dimension, and :math:`d` the number of dimensions.
For example for a two-dimensional quadrilateral element

.. math::

    \underline{f} =
    \big[
        f_x^{(0)},
        f_y^{(0)},
        f_x^{(1)},
        f_y^{(1)},
        f_x^{(2)},
        f_y^{(2)},
        f_x^{(3)},
        f_y^{(3)}
    \big]^T

