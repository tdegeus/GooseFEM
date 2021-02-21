
*********
Changelog
*********

v0.9.0
======

API Changes
-----------

*   VectorPartitioned::asDofs(dofval_u, dofval_p, dofval) ->
    VectorPartitioned::dofsFromParitioned(dofval_u, dofval_p, dofval)

*   VectorPartitioned::asNode(dofval_u, dofval_p, nodevec) ->
    VectorPartitioned::nodeFromPartitioned(dofval_u, dofval_p, nodevec)

*   VectorPartitioned::asElement(dofval_u, dofval_p, elemvec) ->
    VectorPartitioned::elementFromPartitioned(dofval_u, dofval_p, elemvec)

*   Version defines as replaced by ``#define GOOSEFEM_VERSION``,
    added convenience methods ``GooseFEM::version()`` and ``GooseFEM::version_dependencies()``.

Deprecating in next version
----------------------------

*   VectorPartitioned::assembleDofs_u
*   VectorPartitioned::assembleDofs_p
*   Mesh::Renumber::get
*   Mesh::Reordered::get

Changes under the hood
----------------------

*   Overloading from Vector (also in Python API)
*   Overloading from QuadratureBase (also in Python API)
*   Added doxygen docs (published to GitHub pages)

v0.8.6
======

*   String-define safety: stringification + unquoting.

v0.8.2
======

*   Using setuptools_scm to manage version (#169)

v0.8.1
======

*   Various documentation updates: using doxygen (e.g. #168, #167, #157, #150)
    *  Adding autodocs using doxygen/breathe.
    *  Adding autodocs Python API with references to the C++ docs.
*   Using GitHub pages for doxygen docs (#156, #155)
*   Adding version information (incl. git commit hash) (#166)
*   Adding GooseFEM::Element::Quad4::Quadrature::interp_N_vector
*   Generalizing GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular
*   Generalising implementation:
    *   Internally deriving from Vector
    *   Python API: unifying Element
    *   Python API: fixing overloaded methods
    *   Removing internal use of deprecated method
    *   Using "initQuadratureBase" in derived Quadrature classes
    *   Introducing QuadratureBase class -> avoids copies of convenience functions
*   [CI] Using ctest command to improve output in case of test failure
*   Restructuring environment (#154)
*   Fixing readthedocs setup (#153)

v0.8.0
======

*   [CI] Using gcc-8
*   Adding Mesh::Quad4::FineLayer::elementsLayer
*   Stitch: Adding nodesets to example
*   Stitch: Adding hybrid example. Adding assertions.
*   Making API more functional
*   Adding Mesh::ManualStich
*   Adding Mesh::Stitch
*   Minor style update
*   [CMake] Minor updates in testing
*   [CI] improve comments (#142)
*   Combining tests MeshQuad4 (#141)
*   Using clang on Windows (#139)

v0.7.0
======

*   Adding ``Mesh::Quad4::FineLayer::elementgrid_leftright``

v0.6.1
======

*   Minor bugfix ``Mesh::Quad4::FineLayer::elementgrid_around_ravel``: allowing huge sizes.

v0.6.0
======

*   Adding ``Mesh::Quad4::FineLayer::elementgrid_around_ravel``
*   ``FineLayer::elementgrid_ravel``: Adding test
*   Renaming ``elementMatrix`` -> ``elementgrid`` everywhere
*   Adding ``Mesh::Quad4::FineLayer::elementgrid_ravel``
*   Adding ``GOOFEM_WIP_ASSERT`` to assert if code needs to be generalized
*   API change: renaming ``Mesh::Quad4::Regular::elementMatrix``
    -> M``esh::Quad4::Regular::elementgrid``.

v0.5.1
======

*   FineLayer - replica: bug-fix in size detection.
*   Updated examples to new GMat API.

v0.5.0
======

*   Renaming ``MatrixDiagonal::AsDiagonal`` -> ``MatrixDiagonal::Todiagonal``
    to maintain API consistency.
*   Adding ``Mesh::elemmap2nodemap``. Updating Python API.
*   Adding ``roll`` to FineLayer.
*   Adding ``Mesh::centers`` and ``Mesh::defaultElementType``.
*   Mapping connectivity on generating FineLayer-object.
*   Switching to new GMat API.
*   Solver: force factorization on the first call.
*   Sorting output of ``GooseFEM::Mesh::elem2node``. Adding checks.
*   Switched to GitHub CI.
*   Adding ``todense`` to sparse matrix classes.
*   Adding ``dot`` to ``MatrixPartitioned``.

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
