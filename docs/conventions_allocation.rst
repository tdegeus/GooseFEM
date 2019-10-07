.. _conventions_allocation:

Data allocation
===============

Most of GooseFEM's functions provided an interface to:

*   Allocate output arrays: names start with a **upper-case** letter. For example:

    .. code-block:: cpp

        ue = GooseFEM::Vector::AsElement(disp);

*   Directly write to output arrays, without allocation them and copying them, by taking a pointer the externally allocated array as the last input argument(s): names start with a **lower-case** letter. For example:

    .. code-block:: cpp

        GooseFEM::Vector::asElement(disp, ue);


