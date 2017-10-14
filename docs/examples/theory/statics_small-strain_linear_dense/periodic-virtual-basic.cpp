/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <iostream>
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseMaterial/Metal/LinearStrain/Elastic.h>

// alias relevant types for cppmat - exclusively 3-D
using T4  = cppmat::cartesian3d::tensor4 <double>;
using T2  = cppmat::cartesian3d::tensor2 <double>;
using T2s = cppmat::cartesian3d::tensor2s<double>;

// alias relevant types for Eigen
using MatD = GooseFEM::MatD;
using MatS = GooseFEM::MatS;
using ColD = GooseFEM::ColD;
using ColS = GooseFEM::ColS;

// alias
namespace GooseElas = GooseMaterial::Metal::LinearStrain::Elastic;

// =================================================================================================

int main()
{

  // --------
  // geometry
  // --------

  // create a mesh
  GooseFEM::Mesh::Quad4::Regular mesh(3,3);

  // element type for extrapolation and quadrature
  GooseFEM::Quad4 el;

  // extract relevant mesh data
  size_t nnode = mesh.nnode(); // number of nodes
  size_t nelem = mesh.nelem(); // number of elements
  size_t ndim  = mesh.ndim (); // number of dimensions        (== 2)
  size_t nne   = mesh.nne  (); // number of nodes per element (== 4)
  MatS   conn  = mesh.conn (); // connectivity (nodes of each element) : one row per element
  MatD   x0    = mesh.coor (); // nodal position in x- and y-direction : one row per node
  size_t ndof  = nnode * ndim; // total number of DOFs

  // add virtual nodes to prescribe the macroscopic deformation
  // - node numbers
  ColS nodesVirtual ( ndim );
  for ( size_t i = 0 ; i < ndim ; ++i ) nodesVirtual(i) = nnode + i;
  // - extend system
  nnode += ndim;
  ndof  += ndim * ndim;
  // - extend nodal position
  x0.conservativeResize(nnode,Eigen::NoChange);
  // - assign dummy positions (to facilitate further use)
  for ( size_t i = 0 ; i < ndim ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      x0(nodesVirtual(i),j) = 0.0;

  // nodal displacements
  // - allocate as a column of vectors; like "x0"
  MatD u(nnode,ndim);
  // - zero-initialize
  u.setZero();

  // DOF-numbers for each vector component of each node (sequential)
  MatS dofs = GooseFEM::Mesh::dofs(nnode,ndim);

  // -----------
  // periodicity
  // -----------

  // allocate size of the independent and dependent parts of the system and tying matrix
  size_t nni,nnd;
  MatD   C_di;

  // renumber DOFs in the order [ iii , iid ]; create matrix with nodal tyings
  {
    // - get periodic nodes
    MatS nodesPeriodic = mesh.nodesPeriodic();
    // - sizes
    nnd = ndim * nodesPeriodic.rows();
    nni = ndof - nnd;
    // - allocate lists of dependent/independent DOFs
    ColS iid(nnd);
    // - extract dependent DOFs
    for ( size_t i = 0 ; i < nodesPeriodic.rows() ; ++i )
      for ( size_t j = 0 ; j < ndim ; ++j )
        iid(i*ndim+j) = dofs(nodesPeriodic(i,1),j);
    // - renumber DOFs in the following order [ iii , iid ]
    dofs = GooseFEM::Mesh::renumber( dofs , iid , "end" );

    // - allocate matrix with nodal tyings
    C_di.conservativeResize(nnd,nni);
    // - zero initialize
    C_di.setZero();
    // - fill
    for ( size_t i = 0 ; i < nodesPeriodic.rows() ; ++i )
    {
      // -- get dependent ("dep") and independent ("ind") node numbers
      size_t nd_dep = nodesPeriodic(i,1);
      size_t nd_ind = nodesPeriodic(i,0);
      // -- tie all DOFs of those nodes
      for ( size_t j = 0 ; j < ndim ; ++j )
      {
        // --- get relevant DOF numbers
        size_t df_dep = dofs(nd_dep,j);
        size_t df_ind = dofs(nd_ind,j);
        // --- u_dep = u_ind + ...
        C_di(  df_dep-nni , df_ind ) = +1.0;
        // --- u_dep = ... + (F-I).dot(X_dep - X_ind)
        for ( size_t k = 0 ; k < ndim ; ++k )
          C_di( df_dep-nni , dofs(nodesVirtual(j),k) ) = x0(nd_dep,k) - x0(nd_ind,k);
      }
    }
  }

  // -------------------
  // boundary conditions
  // -------------------

  // allocate size of the partitioned parts
  size_t nnp = ( ndim+1 ) * ndim;
  size_t nnu = nni - nnp;
  // - allocate lists of prescribed/unknown DOFs and their values
  ColS iip (nnp);
  ColS iiu (nnu);
  ColD DU_p(nnp);

  // set prescribed/unknown DOFs and their values
  {
    // - reference node to suppress rigid body motion
    size_t nodeRef = mesh.nodesRef();

    // - macroscopic deformation gradient update
    // -- allocate
    T2 Fbar(0.0);
    // -- set non-zeros values
    Fbar(0,1) += 0.01;

    // - prescribe displacement-update DOFs and their values
    // -- counter
    size_t k = 0;
    // -- macroscopic deformation
    for ( size_t i = 0 ; i < ndim ; ++i ) {
      for ( size_t j = 0 ; j < ndim ; ++j ) {
        iip(k) = dofs(nodesVirtual(i),j); DU_p(k) = Fbar(i,j); ++k;
      }
    }
    // -- suppress rigid body motion
    for ( size_t j = 0 ; j < ndim ; ++j ) {
      iip(k) = dofs(nodeRef,j); DU_p(k) = 0.; ++k;
    }
    // -- list with all 'independent' DOFs
    ColS dofs_p = ColS::LinSpaced(nni, 0, nni);
    // -- temporary copy of "iip"
    ColS tmp = iip;
    // -- sort "iip"
    std::sort ( tmp.data(), tmp.data()+nnp );
    // -- "iiu" = setdiff ( "dofs_p" , "iip" )
    std::set_difference ( dofs_p.data(),dofs_p.data()+nni,tmp.data(),tmp.data()+nnp,iiu.data() );
  }

  // --------
  // material
  // --------

  std::vector<GooseElas::Cartesian3d::Material> mat;

  for ( size_t i = 0 ; i < 1 ; ++i )
  {
    // set elastic parameters
    double E  = 100.0;
    double nu =   0.3;
    double kappa,G;
    // convert elastic parameters to the pair required by "GooseElas::Cartesian3d::Material"
    std::tie(kappa,G) = GooseElas::ConvertParameters("E,nu",E,nu,"K,G");
    // element material definition: linear elasticity
    GooseElas::Cartesian3d::Material emat = GooseElas::Cartesian3d::Material(kappa,G);
    // add to list
    mat.push_back(emat);
  }

  for ( size_t i = 1 ; i < nelem ; ++i )
  {
    // set elastic parameters
    double E  =   1.0;
    double nu =   0.3;
    double kappa,G;
    // convert elastic parameters to the pair required by "GooseElas::Cartesian3d::Material"
    std::tie(kappa,G) = GooseElas::ConvertParameters("E,nu",E,nu,"K,G");
    // element material definition: linear elasticity
    GooseElas::Cartesian3d::Material emat = GooseElas::Cartesian3d::Material(kappa,G);
    // add to list
    mat.push_back(emat);
  }

  // -------------------------
  // global system - variables
  // -------------------------

  // global system - partitioned in independent and dependent DOFs
  MatD K_ii(nni,nni); // stiffness           (matrix)
  MatD K_id(nni,nnd); // "
  MatD K_di(nnd,nni); // "
  MatD K_dd(nnd,nnd); // "
  ColD F_i (nni    ); // internal force      (column)
  ColD F_d (nnd    ); // "
  ColD DU_i(nni    ); // displacement update (column)
  ColD DU_d(nnd    ); // "

  // --------------------------------------
  // element & quadrature-point - variables
  // --------------------------------------

  // quadrature-point tensors
  // - displacement gradient (2-D)
  MatD gradu(ndim,ndim);
  // - strain, stress, and stiffness (3-D; strain has zero z-components)
  T2s eps(0.0), sig;
  T4  K4;

  // element arrays
  MatD k_e  (nne*ndim,nne*ndim); // element stiffness (matrix)
  ColD f_e  (nne*ndim         ); // element internal force (column)
  ColS dof_e(nne*ndim         ); // array with DOF numbers of each row/column of "f_e" and "k_e"
  MatD x0_e (nne     ,    ndim); // nodal positions     (column of vectors)
  MatD u_e  (nne     ,    ndim); // nodal displacements (column of vectors)

  // ------------------------
  // global system - assemble
  // ------------------------

  // zero-initialize
  K_ii.setZero();
  K_id.setZero();
  K_di.setZero();
  K_dd.setZero();
  F_i .setZero();
  F_d .setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers; nodal positions and displacements for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) u_e  .row(m)               = u   .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - zero initialize element stiffness and internal force
    k_e.setZero();
    f_e.setZero();

    // - loop over the quadrature-points
    for ( size_t k = 0 ; k < el.QuadGaussNumPoints() ; ++k )
    {
      // -- evaluate the (gradient of the) shape functions
      el.eval( x0_e, k );

      // -- local displacement gradient
      gradu = el.dNdx.transpose() * u_e;

      // -- local strain tensor (stored symmetrically, 2-D plane strain: all z-components zero)
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = i ; j < ndim ; ++j )
          eps(i,j) = 0.5 * ( gradu(i,j) + gradu(j,i) );

      // -- local constitutive response: stiffness and stress tensors
      std::tie(K4,sig) = mat[e].tangent_stress(eps);

      // -- add stress and stiffness tensors to element internal force column and stiffness matrix
      f_e += el.gradN_tensor2(sig);
      k_e += el.gradN_tensor4_gradNT(K4);
    }

    // - assemble element internal force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) {
      if ( dof_e(i) < nni ) F_i ( dof_e(i)     ) += f_e ( i );
      else                  F_d ( dof_e(i)-nni ) += f_e ( i );
    }

    // - assemble element stiffness to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) {
      for ( size_t j = 0 ; j < nne*ndim ; ++j ) {
        if      ( dof_e(i) < nni and dof_e(j) < nni ) K_ii(dof_e(i)    , dof_e(j)    ) += k_e(i,j);
        else if ( dof_e(i) < nni                    ) K_id(dof_e(i)    , dof_e(j)-nni) += k_e(i,j);
        else if (                    dof_e(j) < nni ) K_di(dof_e(i)-nni, dof_e(j)    ) += k_e(i,j);
        else                                          K_dd(dof_e(i)-nni, dof_e(j)-nni) += k_e(i,j);
      }
    }
  }

  // ------------------------
  // eliminate dependent DOFs
  // ------------------------

  // eliminate dependent DOFs from the system
  K_ii += K_id * C_di + C_di.transpose() * K_di + C_di.transpose() * K_dd * C_di;
  F_i  +=               C_di.transpose() * F_d;

  // -------------------
  // partition and solve
  // -------------------

  // allocate partitioned system
  MatD K_uu(nnu,nnu);
  MatD K_up(nnu,nnp);
  MatD K_pu(nnp,nnu);
  MatD K_pp(nnp,nnp);
  ColD F_u (nnu    );
  ColD DU_u(nnu    );

  // partition system
  for ( size_t i = 0 ; i < nnu ; ++i )
    for ( size_t j = 0 ; j < nnu ; ++j )
      K_uu(i,j) = K_ii( iiu(i) , iiu(j) );
  // -
  for ( size_t i = 0 ; i < nnu ; ++i )
    for ( size_t j = 0 ; j < nnp ; ++j )
      K_up(i,j) = K_ii( iiu(i) , iip(j) );
  // -
  for ( size_t i = 0 ; i < nnp ; ++i )
    for ( size_t j = 0 ; j < nnu ; ++j )
      K_pu(i,j) = K_ii( iip(i) , iiu(j) );
  // -
  for ( size_t i = 0 ; i < nnp ; ++i )
    for ( size_t j = 0 ; j < nnp ; ++j )
      K_pp(i,j) = K_ii( iip(i) , iip(j) );
  // -
  for ( size_t i = 0 ; i < nnu ; ++i ) F_u ( i ) = F_i ( iiu(i) );

  // solve for the unknown DOFs
  DU_u = K_uu.ldlt().solve( - F_u - K_up * DU_p );

  // place to global system of independent DOFs
  for ( size_t i = 0 ; i < nnu ; ++i ) DU_i( iiu(i) ) = DU_u(i);
  for ( size_t i = 0 ; i < nnp ; ++i ) DU_i( iip(i) ) = DU_p(i);

  // reconstruct the displacement of the dependent DOFs
  DU_d = C_di * DU_i;

  // add displacement update to global system
  for ( size_t i = 0 ; i < ndof ; ++i ) {
    if ( dofs(i) < nni ) u ( i ) += DU_i ( dofs(i)       );
    else                 u ( i ) += DU_d ( dofs(i) - nni );
  }

  // print the result to the screen
  std::cout << "new nodal coordinates = " << std::endl << x0 + u << std::endl << std::endl;

  // ----------------------
  // compute reaction force
  // ----------------------

  // zero-initialize
  F_i.setZero();
  F_d.setZero();

  // also compute volume averaged stress tensor, to show that is the same as the reaction force
  T2 sigbar(0.0);

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers, nodal positions, and displacements for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) u_e  .row(m)               = u   .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - zero initialize element and internal force
    f_e.setZero();

    // - loop over the quadrature-points
    for ( size_t k = 0 ; k < el.QuadGaussNumPoints() ; ++k )
    {
      // -- evaluate the (gradient of the) shape functions
      el.eval( x0_e, k );

      // -- local displacement gradient
      gradu = el.dNdx.transpose() * u_e;

      // -- local strain tensor (stored symmetrically, 2-D plane strain: all z-components zero)
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = i ; j < ndim ; ++j )
          eps(i,j) = 0.5 * ( gradu(i,j) + gradu(j,i) );

      // -- local constitutive response: stress tensor
      sig = mat[e].stress(eps);

      // -- add stress tensor to element internal force column
      f_e += el.gradN_tensor2(sig);

      // -- add to volume averaged stress tensor
      sigbar += sig * el.V;
    }

    // - assemble internal force column to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) {
      if ( dof_e(i) < nni ) F_i( dof_e(i)     ) += f_e(i);
      else                  F_d( dof_e(i)-nni ) += f_e(i);
    }
  }

  // eliminate dependent DOFs
  F_i += C_di.transpose() * F_d;

  // show macroscopic reaction force
  MatD Tbar(ndim,ndim);

  for ( size_t i = 0 ; i < ndim ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      Tbar(i,j) = F_i ( dofs(nodesVirtual(i),j) );

  // print to screen
  std::cout << "Tbar   = " << std::endl << Tbar << std::endl;
  std::cout << "sigbar = " << std::endl;
  sigbar.printf("%16.8f");

  return 0;
}
