
*************
Semi-periodic
*************

Consider a "rheometer" which corresponds to a torus, in which a linear elastic material is placed
consists of two phases.
In terms of boundary conditions this implies taking the geometry periodic in horizontal direction,
while the displacements of the top and bottom boundary are controlled.
For simplicity two-dimensional plane strain is considered.

| :download:`sketch.svg <statics/partial-periodic_elastic_sketch.svg>`

Below the important difference with respect to the previous example are discussed.
The full example can be downloaded:

:download:`example.py <statics/partial-periodic_elastic.py>`

Apply periodicity
=================

.. literalinclude:: statics/partial-periodic_elastic.py
    :language: py
    :lines: 27-28
    :emphasize-lines: 1

Applying periodicity in this case is rather straightforward.
In particular the degrees-of-freedom along the right edge are eliminated,
and replaced by the degrees-of-freedom of the left edge.
The size of the actually solved system is therefore reduced, while the response vectors are
simply assembled to both sides of the geometry.

Fixed displacement
==================

.. literalinclude:: statics/partial-periodic_elastic.py
    :language: py
    :lines: 86

The degrees-of-freedom of which the displacement is controlled are finally extracted
from the renumbered list of degrees-of-freedom.
