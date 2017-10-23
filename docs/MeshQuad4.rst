
**********************
<GooseFEM/MeshQuad4.h>
**********************

``GooseFEM::Mesh::Quad4::Regular``
==================================

Regular mesh of linear quadrilaterals in two-dimensions.

Methods:

*  ``MatS = GooseFEM::Mesh::Quad4::Regular.nodesPeriodic()``

    Returns a list of periodic node-pairs in which the first column contains the independent nodes, and the second column the dependent nodes.


``GooseFEM::Mesh::Quad4::FineLayer``
====================================

Regular mesh with a fine layer of quadrilateral elements, and coarser elements above and below.

.. image:: figures/MeshQuad4/FineLayer/example.svg
  :width: 500px
  :align: center

.. note::

  The coarsening depends strongly on the desired number of elements in horizontal elements. The becomes clear from the following example:

  .. code-block:: python

    mesh = gf.Mesh.Quad4.FineLayer(6*9  ,51) # left   image :  546 elements
    mesh = gf.Mesh.Quad4.FineLayer(6*9+3,51) # middle image :  703 elements
    mesh = gf.Mesh.Quad4.FineLayer(6*9+1,51) # right  image : 2915 elements

  .. image:: figures/MeshQuad4/FineLayer/behavior.svg
    :width: 1000px
    :align: center

Methods:

*   ``MatS = GooseFEM::Mesh::Quad4::Regular.nodesPeriodic()``

    Returns a list of periodic node-pairs in which the first column contains the independent nodes, and the second column the dependent nodes.

*   ``ColS = GooseFEM::Mesh::Quad4::Regular.elementsFine()``

    Returns a list with the element numbers of the fine elements in the center of the mesh (highlighted in the plot below).

    .. image:: figures/MeshQuad4/FineLayer/example_elementsFine.svg
      :width: 500px
      :align: center
