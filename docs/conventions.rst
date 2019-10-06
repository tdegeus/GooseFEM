
***********
Conventions
***********

Data-storage
============

+-----------+------------------------------------------------+----------------------------------+---------------------------+
|  Alias    | Description                                    | Shape                            | Type                      |
+===========+================================================+==================================+===========================+
| "dofval"  | degrees-of-freedom                             | [ndof]                           | ``xt::xtensor<double,1>`` |
+-----------+------------------------------------------------+----------------------------------+---------------------------+
| "nodevec" | nodal vectors                                  | [nnode, ndim]                    | ``xt::xtensor<double,2>`` |
+-----------+------------------------------------------------+----------------------------------+---------------------------+
| "elemvec" | nodal vectors stored per element               | [nelem, nne, ndim]               | ``xt::xtensor<double,3>`` |
+-----------+------------------------------------------------+----------------------------------+---------------------------+
| "elemmat" | matrices stored per element                    | [nelem, nne*ndim, nne*ndim]      | ``xt::xtensor<double,3>`` |
+-----------+------------------------------------------------+----------------------------------+---------------------------+
| "qtensor" | tensors stored (as list) per integration point | [nelem, nip, #tensor-components] | ``xt::xtensor<double,4>`` |
+-----------+------------------------------------------------+----------------------------------+---------------------------+
| "qscalar" | scalars stored per integration point           | [nelem, nip]                     | ``xt::xtensor<double,2>`` |
+-----------+------------------------------------------------+----------------------------------+---------------------------+

Data-allocation
===============

Most of GooseFEM's functions provided an interface to:

*   Allocate output arrays: names start with a **upper-case** letter. For example:

    .. code-block:: cpp

        ue = GooseFEM::Vector::AsElement(disp);

*   Directly write to output arrays, without allocation them and copying them, by taking a pointer the externally allocated array as the last input argument(s): names start with a **lower-case** letter. For example:

    .. code-block:: cpp

        GooseFEM::Vector::asElement(disp, ue);


