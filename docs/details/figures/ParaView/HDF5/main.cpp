#include <GooseFEM/GooseFEM.h>
#include <GooseFEM/ParaView.h>

namespace PV = GooseFEM::ParaView::HDF5;

int main()
{
    // define mesh
    GooseFEM::Mesh::Quad4::FineLayer mesh(6, 18);

    // extract mesh fields
    xt::xtensor<double, 2> coor = mesh.coor();
    xt::xtensor<double, 2> conn = mesh.conn();
    xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());

    // vector definition:
    // provides methods to switch between dofval/nodeval/elemvec, or to manipulate a part of them
    GooseFEM::Vector vector(conn, mesh.dofs());

    // FEM quadrature
    GooseFEM::Element::Quad4::Quadrature elem(vector.AsElement(coor));

    // open output file
    H5Easy::File data("output.h5", H5Easy::File::Overwrite);

    // initialise ParaView metadata
    PV::TimeSeries xdmf;

    // save mesh to output file
    H5Easy::dump(data, "/coor", coor);
    H5Easy::dump(data, "/conn", conn);

    // define strain history
    xt::xtensor<double, 1> gamma = xt::linspace<double>(0, 1, 100);

    // loop over increments
    for (size_t inc = 0; inc < gamma.size(); ++inc) {
        // apply fictitious strain
        for (size_t node = 0; node < disp.shape(0); ++node)
            disp(node, 0) = gamma(inc) * (coor(node, 1) - coor(0, 1));

        // compute strain tensor
        xt::xtensor<double, 4> Eps = elem.SymGradN_vector(vector.AsElement(disp));
        xt::xtensor<double, 1> Eps_xy = xt::view(Eps, xt::all(), 0, 0, 1);

        // store data to output file
        H5Easy::dump(data, "/disp/" + std::to_string(inc), PV::as3d(disp));
        H5Easy::dump(data, "/eps_xy/" + std::to_string(inc), Eps_xy);

        // add increment to ParaView metadata
        xdmf.push_back(PV::Increment(
            PV::Connectivity(data, "/conn", mesh.getElementType()),
            PV::Coordinates(data, "/coor"),
            {PV::Attribute(
                 data, "/disp/" + std::to_string(inc), "Displacement", PV::AttributeType::Node),
             PV::Attribute(
                 data, "/eps_xy/" + std::to_string(inc), "Eps_xy", PV::AttributeType::Cell)}));
    }

    // write ParaView metadata
    xdmf.write("output.xdmf");

    return 0;
}
