
*********
Changelog
*********

v0.5.0
======

*   Renaming ``MatrixDiagonal::AsDiagonal`` -> ``MatrixDiagonal::Todiagonal``
    to maintain API consistency.
*   Adding ``Mesh::elemmap2nodemap``. Updating Python API.
*   Adding ``roll`` to FineLayer
*   Adding ``Mesh::centers`` and ``Mesh::defaultElementType``
*   Mapping connectivity on generating FineLayer-object
*   Switching to new GMat API

v0.4.2
======

*   CMake: using Eigen's CMake target.

v0.4.1
======

API additions
-------------

*   Added  "AllocateElemmat".

v0.4.0
======

API additions
-------------

*   Added "AllocateQtensor", "AllocateQscalar", "AllocateDofval", "AllocateNodevec", "AllocateElemvec".

API changes
-----------

*   Removing Paraview interface: replaced by external libraries "XDMFWrite_HighFive" and "XDMFWrite_h5py".

*   Element*: "dV" now only returns raw data, the "asTensor" member function (and free function) can be used to convert the 'qscalar' to a 'qtensor'.

*   Separating sparse solver in separate class to offer more flexibility in the future.

*   Adding "dot" to "Matrix".

Other updates
-------------

*   Applying clang-format to source, python API, tests, and examples..

*   Adding test GMatElastoPlasticQPot.

*   Adding test based on hybrid material definitions.

*   Formatting update: renaming all return variables "out" to "ret".

*   Correction zero allocation to allows for dofval.size() > nodevec.size()

*   Formatting update xt::amax and xt::sum.

*   Renaming private function to begin with caps when the function allocates its return data.

*   Reducing copies when using Eigen.

*   Reducing default size examples.

*   Supporting Windows (#87).

*   Removing xtensor_fixed.

*   Using xt::has_shape.
