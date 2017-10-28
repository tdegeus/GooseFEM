#include <GooseFEM/GooseFEM.h>

using ColS = GooseFEM::ColS;
using MatS = GooseFEM::MatS;

int main()
{
  // mesh definition
  GooseFEM::Mesh::Quad4::Regular mesh = GooseFEM::Mesh::Quad4::Regular(2,2);

  // extract information
  size_t ndim   = mesh.ndim       ();
  ColS   top    = mesh.nodesTop   ();
  ColS   bottom = mesh.nodesBottom();
  ColS   left   = mesh.nodesLeft  ();
  ColS   right  = mesh.nodesRight ();
  size_t nleft  = static_cast<size_t>(left.size());
  size_t ntop   = static_cast<size_t>(top .size());

  // default DOF-numbers: sequential numbering per element
  MatS dofs = mesh.dofs();

  std::cout << "dofs = " << std::endl;
  std::cout << dofs << std::endl << std::endl;

  // periodicity in horizontal direction : eliminate 'dependent' DOFs
  // N.B. the corner nodes are part of the fixed boundaries
  for ( size_t i = 1 ; i < nleft-1 ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      dofs(right(i),j) = dofs(left(i),j);

  // renumber "dofs" to be sequential
  dofs = GooseFEM::Mesh::renumber(dofs);

  std::cout << "periodicity :" << std::endl << "dofs = " << std::endl;
  std::cout << dofs << std::endl << std::endl;

  // construct list with prescribed DOFs
  // - allocate
  ColS fixedDofs(2*ntop*ndim);
  // - read: bottom
  for ( size_t i = 0 ; i < ntop ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      fixedDofs(i*ndim+j) = dofs(bottom(i),j);
  // - read: top
  for ( size_t i = 0 ; i < ntop ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      fixedDofs(i*ndim+j+ntop*ndim) = dofs(top(i),j);

  // renumber "dofs" to have the "fixedDofs" at the end
  dofs = GooseFEM::Mesh::renumber(dofs,fixedDofs);

  std::cout << "fixedDofs to end :" << std::endl << "dofs = " << std::endl;
  std::cout << dofs << std::endl << std::endl;

}
