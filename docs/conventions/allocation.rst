Data allocation
===============

Most of GooseFEM's functions provided an interface to:

1.  Allocate output arrays: names start with a **upper-case** letter. For example:

    .. code-block:: cpp

        ue = GooseFEM::Vector::AsElement(disp);

2.  Directly write to output arrays, without allocation them and copying them, by taking a pointer the externally allocated array as the last input argument(s): names start with a **lower-case** letter. For example:

    .. code-block:: cpp

        GooseFEM::Vector::asElement(disp, ue);

.. note::

    The Python API only provides option 1. Option 2 is only available in the C++ API.
