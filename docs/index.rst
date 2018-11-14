
*********
GooseFEM
*********

.. |badge1| image:: https://img.shields.io/badge/license-GPLv3-brightgreen.svg
  :target: https://github.com/tdegeus/GooseFEM/blob/master/LICENSE
  :alt: GPLv3 license

.. |badge2| image:: https://img.shields.io/badge/warranty-no-red.svg
  :target: https://github.com/tdegeus/GooseFEM/blob/master/LICENSE
  :alt: GPLv3 license

.. |badge3| image:: https://img.shields.io/badge/download-.zip-lightgray.svg
  :target: https://github.com/tdegeus/GooseFEM/zipball/master
  :alt: Download as .zip

.. |badge4| image:: https://img.shields.io/badge/download-.tar.gz-lightgray.svg
  :target: https://github.com/tdegeus/GooseFEM/tarball/master
  :alt: Download as .tar.gz

.. |badge5| image:: https://img.shields.io/badge/contact-tom@geus.me-blue.svg
  :target: mailto:tom@geus.me
  :alt: Contact tom@geus.me

.. |badge6| image:: https://img.shields.io/badge/contact-www.geus.me-blue.svg
  :target: http://www.geus.me
  :alt: Website www.geus.me

.. |badge7| image:: https://img.shields.io/badge/GitHub-tdegeus/GooseFEM-blue.svg
  :target: https://github.com/tdegeus/GooseFEM
  :alt: Github tdegeus/GooseFEM

.. |badge8| image:: https://img.shields.io/badge/documentation-GooseFEM.geus.me-blue.svg
  :target: http://GooseFEM.geus.me
  :alt: Website GooseFEM.geus.me

| |badge1| |badge2| |badge3| |badge4|
| |badge5| |badge6| |badge7|
| |badge8|

.. note::

  This library is free to use under the `GPLv3 license <https://github.com/tdegeus/GooseFEM/blob/master/LICENSE>`_. Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word of mouth advertisement, .... Bugs or feature requests can be filed on `GitHub <http://github.com/tdegeus/GooseFEM>`_. As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.

.. tip::

  This document should be considered as a quick-start guide. A lot effort has been spent on the readability of the code itself (in particular the ``*.h`` files should be instructive). One is highly encouraged to answer more advanced questions that arise from this guide directly using the code. Download buttons to the relevant files are included throughout this reader.

This header-only module provides C++ classes and several accompanying methods to work with n-d arrays and/or tensors. It's usage, programmatically and from a compilation perspective, is really simple. One just has to ``#include <GooseFEM/GooseFEM.h>`` and tell your compiler where GooseFEM is located (and to use the C++14 or younger standard). Really, that's it!

Data-types
==========

[:download:`GooseFEM/GooseFEM.h <../include/GooseFEM/GooseFEM.h>`]

Beyond the default C++ types, GooseFEM use Eigen columns/matrices and cppmat multi-dimensional arrays. In particular:

+------------+---------------------------------------------+
| Type       | Description                                 |
+============+=============================================+
| ``MatD``   | row-major ``Eigen::matrix``                 |
+------------+---------------------------------------------+
| ``ColD``   | (column-major) ``Eigen::Matrix`` (column)   |
+------------+---------------------------------------------+
| ``SpMatD`` | row-major ``Eigen::SparseMatrix``           |
+------------+---------------------------------------------+
| ``ArrD``   | multi-dimensional ``cppmat::array``         |
+------------+---------------------------------------------+

The last letter thereby indicates the type specialisation:

+------------+----------------+
| Indicator  | Description    |
+============+================+
| ``D``      | ``double``     |
+------------+----------------+
| ``S``      | ``size_t``     |
+------------+----------------+
| ``I``      | ``int``        |
+------------+----------------+

Data-storage
============

+-----------+------------------------------------------------+-------------+----------------------------------+
|  Alias    | Description                                    | Type        | Shape                            |
+===========+================================================+=============+==================================+
| "dofval"  | degrees-of-freedom                             | ``ColD``    | [ndof]                           |
+-----------+------------------------------------------------+-------------+----------------------------------+
| "nodevec" | nodal vectors                                  | ``MatD``    | [nnode, ndim]                    |
+-----------+------------------------------------------------+-------------+----------------------------------+
| "elemvec" | nodal vectors stored per element               | ``ArrD``    | [nelem, nne, ndim]               |
+-----------+------------------------------------------------+-------------+----------------------------------+
| "elemmat" | matrices stored per element                    | ``ArrD``    | [nelem, nne*ndim, nne*ndim]      |
+-----------+------------------------------------------------+-------------+----------------------------------+
| "qtensor" | tensors stored (as list) per integration point | ``ArrD``    | [nelem, nip, #tensor-components] |
+-----------+------------------------------------------------+-------------+----------------------------------+
| "qscalar" | scalars stored per integration point           | ``ArrD``    | [nelem, nip]                     |
+-----------+------------------------------------------------+-------------+----------------------------------+

Contents
========

.. toctree::
   :caption: USAGE
   :maxdepth: 1

   Mesh.rst
   Element.rst
   Vector.rst

.. toctree::
   :caption: INSTALLATION
   :maxdepth: 1

   compile.rst
   develop.rst

.. toctree::
   :caption: EXAMPLES
   :maxdepth: 1

   examples/Dynamics/readme.rst

.. tip::

  A compact reader covering the basic theory is available `here <https://github.com/tdegeus/GooseFEM/docs/theory/readme.pdf>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

