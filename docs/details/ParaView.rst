
********
ParaView
********

.. seealso::

  | :download:`GooseFEM/ParaView.h <../../include/GooseFEM/ParaView.h>`
  | :download:`GooseFEM/ParaView.hpp <../../include/GooseFEM/ParaView.hpp>`

.. note::

  This header relies on HDF5 and HighFive as dependencies. If you wish to use ParaView support without making use of HDF5 and HighFive, you have to  define ``GOOSEFEM_NO_HIGHFIVE`` before including ``ParaView.h`` for the first time:

  .. code-block:: cpp

    #define GOOSEFEM_NO_HIGHFIVE
    #include <GooseFEM/ParaView.h>

  In this case the library does not automatically read the shapes of the datasets. Instead you'll have to provide them as ``std::vector<size_t>``.

HDF5
----

TimeSeries
^^^^^^^^^^

A TimeSeries is constructed from a number of Increments. The consecutive Increments are added to the TimeSeries using push_back. The order in which this is done will define the order of the TimeSeries. Each Increment is constructed from a mesh (using Connectivity and Coordinates) and a number of nodal or cell Attributes.

Consider this example

:download:`figures/ParaView/HDF5/main.cpp <figures/ParaView/HDF5/main.cpp>`

.. literalinclude:: figures/ParaView/HDF5/main.cpp
   :language: cpp

.. tip::

  A displacement vector in must be always 3-d in ParaView, even when the mesh is in 2-d. Use the ``GooseFEM::ParaView::HDF5::as3d(...)`` function to convert a matrix of 2-d displacements to a matrix of 3-d displacements.

.. note::

  The Python interface avoids the HDF5 and HighFive dependencies. One therefore has to provide the datasets' shapes. Consider the following Python example:

  :download:`figures/ParaView/HDF5/main.py <figures/ParaView/HDF5/main.py>`

  .. literalinclude:: figures/ParaView/HDF5/main.py
     :language: python
