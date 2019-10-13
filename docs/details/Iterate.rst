.. _Iterate:

*******
Iterate
*******

| :download:`GooseFEM/Iterate.h <../../include/GooseFEM/Iterate.h>`
| :download:`GooseFEM/Iterate.hpp <../../include/GooseFEM/Iterate.hpp>`

Iterate::StopList
=================

Convergence check to check that a residual is smaller than a certain value of a certain number of consecutive steps. The class is constructed with the number of consecutive to enforce the residual.

Iterate::StopList::reset(...)
-----------------------------

Reset all residuals to infinity.

Iterate::StopList::stop(...)
----------------------------

Append the list with residuals, return "true" if all residuals are below the tolerance (argument).
