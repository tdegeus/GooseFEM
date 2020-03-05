************
Introduction
************

Overview
--------

This header-only module provides C++ classes (including Python wrappers)
to write simulations using the Finite Element Method.

.. note::

    This library is free to use under the
    `GPLv3 license <https://github.com/tdegeus/GooseFEM/blob/master/LICENSE>`_.
    Any additions are very much appreciated, in terms of suggested functionality, code,
    documentation, testimonials, word of mouth advertisement, ....
    Bug reports or feature requests can be filed on
    `GitHub <http://github.com/tdegeus/GooseFEM>`_.
    As always, the code comes with no guarantee.
    None of the developers can be held responsible for possible mistakes.

Hallmark feature
----------------

The hallmark feature of GooseFEM is that data is not stored in GooseFEM's classes
but only in multi-dimensional arrays with certain
:ref:`storage conventions <conventions_storage>`.
Consequently one is entirely free to mix and match routines,
from GooseFEM or even from somewhere else.
Consequently the readability of the end-user's code and of the library remains high:
Only the main function combines different ingredients,
there is no interdependence between GooseFEM's classes.

Documentation
-------------

This document should be considered as a quick-start guide.
A lot effort has been spent on the readability of the code itself
(in particular the ``.h`` files should be instructive).
One is highly encouraged to answer more advanced questions that arise from this guide
directly using the code.
Download buttons to the relevant files are included throughout this reader.

A compact reader covering the basic theory is available as :download:`PDF <theory/readme.pdf>`.
