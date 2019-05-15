import h5py
import numpy                  as np
import GooseFEM               as gf
import GooseFEM.ParaView.HDF5 as pv

# define mesh
mesh = gf.Mesh.Quad4.FineLayer(6, 18)

# extract mesh fields
coor = mesh.coor();
conn = mesh.conn();
disp = np.zeros(coor.shape)

# vector definition:
# provides methods to switch between dofval/nodeval/elemvec, or to manipulate a part of them
vector = gf.Vector(conn, mesh.dofs())

# FEM quadrature
elem = gf.Element.Quad4.Quadrature(vector.AsElement(coor))

# open output file
data = h5py.File("output.h5", "w")

# initialise ParaView metadata
xdmf = pv.TimeSeries()

# save mesh to output file
data["/coor"] = coor
data["/conn"] = conn

# define strain history
gamma = np.linspace(0, 1, 100);

# loop over increments
for inc in range(len(gamma)):

  # apply fictitious strain
  for node in range(disp.shape[0]):
    disp[node,0] = gamma[inc] * (coor[node,1] - coor[0,1])

  # compute strain tensor
  Eps = elem.SymGradN_vector(vector.AsElement(disp));
  Eps_xy = Eps[:, 0, 0, 1]

  # store data to output file
  data["/disp/"   + str(inc)] = pv.as3d(disp)
  data["/eps_xy/" + str(inc)] = Eps_xy

  # ParaView metadata
  # - initialise Increment
  xdmf_inc = pv.Increment(
    pv.Connectivity(data.filename, "/conn", pv.ElementType.Quadrilateral, conn.shape),
    pv.Coordinates (data.filename, "/coor"                              , coor.shape),
  )
  # - add attributes to Increment
  dataset = "/disp/" + str(inc)
  xdmf_inc.push_back(pv.Attribute(
    data.filename, dataset, "Displacement", pv.AttributeType.Node, data[dataset].shape))
  # - add attributes to Increment
  dataset = "/eps_xy/" + str(inc)
  xdmf_inc.push_back(pv.Attribute(
    data.filename, dataset, "Eps_xy", pv.AttributeType.Cell, data[dataset].shape))
  # - add Increment to TimeSeries
  xdmf.push_back(xdmf_inc)

# write ParaView metadata
xdmf.write("output.xdmf");
