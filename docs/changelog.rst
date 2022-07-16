
*********
Changelog
*********

v1.2.0
======

API changes
-----------

*   [C++] [Python] Several functions now return a pointer to internal data rather than a copy.
    In the Python-API this changes these functions (e.g. `.foo()`) to properties (e.g. `.foo`).

*   [MeshQuad4.h] Renaming:
    -   ``getCoarseMesh()`` to ``coarseMesh()``
    -   ``getFineMesh()`` to ``fineMesh()``
    -   ``getRegularMesh()`` to ``regularMesh()``
    -   ``getFineLayerMesh()`` to ``fineLayerMesh()``
    -   ``getMap()`` to ``map()``
    -   ``getMapFraction()`` to ``mapFraction()``
    (immediate removal in Python, as there are any changes; functions marked as deprecated in C++).

*   [Iterate.h] Renaming:
    -   ``.get()`` to ``.data()``
    (immediate removal in Python, as there are any changes; functions marked as deprecated in C++).

*   [Mesh.h] Renaming:
    -   ``.getElementType()`` to ``.elementType()``
    (immediate removal in Python, as there are any changes; functions marked as deprecated in C++).

*   [Matrix*.h] Adding pointers to ``.data*()``.
    This deprecates ``Todiagonal``
    (immediate removal in Python, as there are any changes; functions marked as deprecated in C++).

*   [Matrix*.h] Adding allocation free API to Python API.

Other changes
-------------

*   Integrating implementation and headers for maintainability.
    The user is expected to get an overview from the doxygen docs.

*   Using CRTP for `Matrix*`.

v1.1.0
======

*   Adding `Mesh::Quad4::FineLayer::elemrow_type` (#238)

v1.0.0
======

Marking that the code is stable enough to be v1.

v0.13.2
=======

*   Adding to Mesh: nodaltyings, node2dof, and periodic overload of elem2node (#237)

v0.13.1
=======

*   Adding "nodal_mass" (#236).

v0.13.0
=======

*   [docs] Updating Python API links (#234)
*   [docs] Switching to breath-apidoc (#234)
*   Adding missing include (#234)
*   Running pre-commit (#234)
*   [docs] Adding center_of_gravity (#234)
*   [docs] Minor updates (#234)
*   Adding: Mesh::center_of_gravity (#234)
*   Adding pre-commit (#233)
*   Fixing CI
*   Fixing various small typos

v0.12.10
========

*   Fixing small typos (#232)
*   Bugfix MatrixPartioned: fixed reordering mistakes "dot", "set", and "add" (tests were not run by mistake) (#232)
*   [Python] Switching to scikit-build (#232)
*   [CMake] Clean-up (#232)
*   Renaming "test" -> "tests" (#232)
*   [CI] Running tests in release mode

v0.12.9
=======

*   Deprecating ``Iterate::StopList::stop*`` in favour of ``Iterate::StopList::insert_roll``, ``Iterate::StopList::descending``, ``Iterate::StopList::all_less``
*   Fixing rare bug in ``Iterate::StopList::stop``
*   Improving assertion messages (#231)

v0.12.8
=======

*   Fixing typo

v0.12.7
=======

*   Removing unused *pyxtensor* dependency

v0.12.6
=======

*   [CMake] Avoiding setuptools_scm dependency if SETUPTOOLS_SCM_PRETEND_VERSION is defined

v0.12.5
=======

*   [Python] Allow CMake customization using ``CMAKE_ARGS``
*   [CMake] Adding option to build examples
*   [docs] Update style. Adding GH build after release

v0.12.2
=======

*   Element: using shape definitions
*   Fixing warnings
*   Removing xtensor-fixed from examples
*   Templating constructor Mesh::Quad4::FineLayer
*   Using CRTP for Element
*   Using CRTP for Mesh
*   Vector: static & dynamic distribution (#213)
*   Templating Mesh (#211)
*   Templating vector (#209)
*   Removing deprecated methods (#210)
*   Allocate: convert to full-template (#208)
*   [Python] Simplifying implementation Mesh
*   [Python] Element registration of base classes
*   [Python] Starting switch pyxtensor -> xtensor-python (#212, #214)
*   [CI] Switching to micromamba (#207)

v0.12.1
=======

*   Integrating Python in CMake (#206)
*   Minor documentation updates

v0.12.0
=======

*   Add "as3d" and "GooseFEM::AsTensor" to Python API
*   Minor documentation updates
*   Adding virtual destructor to silence warning
*   Documenting MatrixDiagonal and MatrixDiagonalPartitioned
*   Documenting MatrixPartitionedTyings
*   Documenting MatrixPartitioned
*   MatrixPartitioned : deriving from Matrix
*   Documentation: minor style updates
*   Documenting Matrix
*   Documentation: minor style updates
*   Documenting VectorPartitionedTyings
*   Documentation: minor style updates
*   Documenting Iterate
*   Mesh::RegularBase* : converting to abstract classes
*   Adding Mesh::Vstack (wraps ManualStich)
*   Adding simple tests Mesh::Hex8
*   Adding version to GH Pages
*   Mesh::Hex8::Regular and Mesh::Hex8::FineLayer : deriving from Base-class
*   Mesh::Hex8::FineLayer renaming private member

v0.11.0
=======

*   Adding ``Int_N_vector_dV``
*   Renaming ``Interp_N_vector`` -> ``InterpQuad_vector``; generalizing dimensions
*   Adding documentation Mesh::Tri3::Regular (by deriving from Mesh::RegularBase2d)
*   Adding documentation Mesh::Quad4::Regular (by deriving from Mesh::RegularBase2d)
*   Adding documentation Mesh::Quad4::FineLayer (by deriving from Mesh::RegularBase2d)
*   Various documentation updates
*   Minor updates CMake/version (#196)
*   Minor updates setup.py
*   Adding convenience functions (Manual)Stitch (#193)
*   [CI] Minor updates
*   [docs] Updating doxygen-awesome
*   [docs] Updating sphinx+breathe build (#189)

v0.10.3
=======

*   Adding convenience functions to Mesh::Stitch
*   Tyings::Periodic : documenting, adding convenience function, adding assertions
*   Tyings::Control : documenting
*   [Python] Reducing dependency requirements
*   [docs] Various minor updates
*   [docs] Doxygen: Using dark theme

v0.10.2
=======

(Future) API Change
-------------------

Note that a new API is introduced, but the old API is kept under deprecation warnings.

*   Renaming: ``Shape...`` -> ``shape_...``
*   Renaming: ``Allocate...`` -> ``allocate_...``

Changes under the hood
----------------------

*   ElementQuad4/ElementQuad4Axisymmetric/ElementQuad4Planar/ElementHex8: deriving from base class.
*   [Python] Reducing compile time deprecation warnings.
*   Using ``GOOSEFEM_WARNING(...)`` for all deprecation messages.
*   Versioning: avoiding un-needed work-around for conda-forge.
*   [CI] Using micromamba.
*   [CI] Adding dynamics example.
*   Adding headers to all files
*   Documentation updates.
*   Moving doxygen configuration to CMake

v0.10.1
=======

*   Python API: Adding Element::Quad4::MidPoint

v0.10.0
=======

API Changes
-----------

Mesh::Quad4::RefineRegular::mapToCoarse -> Mesh::Quad4::RefineRegular::meanToCoarse,
or Mesh::Quad4::RefineRegular::averageToCoarse

Changes under the hood
----------------------

*   Generalisation Mesh::Quad4::RefineRegular::mapToFine
*   Documentation updates
*   Adding explicit deprecation warnings
*   Element::Quad4: moving inv -> detail::inv

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

New functions
-------------

*   VectorPartitioned::dofs_is_u
*   VectorPartitioned::dofs_is_p

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
