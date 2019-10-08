************
Introduction
************

This header-only module provides C++ classes (including Python wrappers) to write simulations using the Finite Element Method.

.. note::

  This library is free to use under the `GPLv3 license <https://github.com/tdegeus/GooseFEM/blob/master/LICENSE>`_. Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word of mouth advertisement, .... Bugs or feature requests can be filed on `GitHub <http://github.com/tdegeus/GooseFEM>`_. As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.

.. tip::

  This document should be considered as a quick-start guide. A lot effort has been spent on the readability of the code itself (in particular the ``*.h`` files should be instructive). One is highly encouraged to answer more advanced questions that arise from this guide directly using the code. Download buttons to the relevant files are included throughout this reader.

.. todo::

  Hallmark feature is that data is always shared using plain n-d arrays.

.. todo::

  VectorPart. and MatrixPart. avoid having to do convert to dofs

.. tip::

  A compact reader covering the basic theory is available `here <https://github.com/tdegeus/GooseFEM/docs/theory/readme.pdf>`_.
