/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_PARAVIEW_HPP
#define GOOSEFEM_PARAVIEW_HPP

#include "ParaView.h"

namespace GooseFEM {
namespace ParaView {
namespace HDF5 {

inline std::string join(const std::vector<std::string>& lines, const std::string& sep)
{
    if (lines.size() == 1) {
        return lines[0];
    }

    std::string ret = "";

    for (auto line : lines) {
        if (ret.size() == 0) {
            ret += line;
            continue;
        }

        if (line[0] == sep[0]) {
            ret += line;
        }
        else if (ret[ret.size() - 1] == sep[0]) {
            ret += line;
        }
        else {
            ret += sep + line;
        }
    }

    return ret;
}

inline std::string indent(size_t n)
{
    std::string ret = "";

    for (size_t i = 0; i < n; ++i) {
        ret += " ";
    }

    return ret;
}

xt::xtensor<double, 2> as3d(const xt::xtensor<double, 2>& data)
{
    GOOSEFEM_ASSERT(data.shape(1) > 0 && data.shape(1) < 4)

    if (data.shape(1) == 3ul) {
        return data;
    }

    xt::xtensor<double, 2> ret = xt::zeros<double>(std::array<size_t, 2>{data.shape(0), 3ul});

    if (data.shape(1) == 2ul) {
        xt::view(ret, xt::all(), xt::keep(0, 1)) = data;
    }

    if (data.shape(1) == 1ul) {
        xt::view(ret, xt::all(), xt::keep(0)) = data;
    }

    return ret;
}

inline ElementType convert(GooseFEM::Mesh::ElementType type)
{
    if (type == GooseFEM::Mesh::ElementType::Tri3) {
        return ElementType::Triangle;
    }

    if (type == GooseFEM::Mesh::ElementType::Quad4) {
        return ElementType::Quadrilateral;
    }

    if (type == GooseFEM::Mesh::ElementType::Hex8) {
        return ElementType::Hexahedron;
    }

    throw std::runtime_error("Unknown GooseFEM::Mesh::ElementType");
}

inline std::string to_string(ElementType type)
{
    if (type == ElementType::Triangle) {
        return "Triangle";
    }

    if (type == ElementType::Quadrilateral) {
        return "Quadrilateral";
    }

    if (type == ElementType::Hexahedron) {
        return "Hexahedron";
    }

    throw std::runtime_error("Unknown GooseFEM::ParaView::HDF5::ElementType");
}

inline std::string to_string(AttributeType type)
{
    if (type == AttributeType::Cell) {
        return "Cell";
    }

    if (type == AttributeType::Node) {
        return "Node";
    }

    throw std::runtime_error("Unknown GooseFEM::ParaView::HDF5::AttributeType");
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Connectivity::Connectivity(
    const H5Easy::File& data, const std::string& dataset, GooseFEM::Mesh::ElementType type)
{
    init(data, dataset, convert(type));
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Connectivity::Connectivity(
    const std::string& fname, const std::string& dataset, GooseFEM::Mesh::ElementType type)
{
    init(fname, dataset, convert(type));
}
#endif

inline Connectivity::Connectivity(
    const std::string& fname,
    const std::string& dataset,
    GooseFEM::Mesh::ElementType type,
    const std::vector<size_t>& shape)
{
    init(fname, dataset, convert(type), shape);
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Connectivity::Connectivity(
    const H5Easy::File& data, const std::string& dataset, ElementType type)
{
    init(data, dataset, type);
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Connectivity::Connectivity(
    const std::string& fname, const std::string& dataset, ElementType type)
{
    init(fname, dataset, type);
}
#endif

inline Connectivity::Connectivity(
    const std::string& fname,
    const std::string& dataset,
    ElementType type,
    const std::vector<size_t>& shape)
{
    init(fname, dataset, type, shape);
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void
Connectivity::init(const H5Easy::File& data, const std::string& dataset, ElementType type)
{
    m_type = type;
    m_dataset = dataset;

    this->readShape(data);
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void
Connectivity::init(const std::string& fname, const std::string& dataset, ElementType type)
{
    m_type = type;
    m_fname = fname;
    m_dataset = dataset;

    H5Easy::File data(m_fname, H5Easy::File::ReadOnly);

    this->readShape(data);
}
#endif

inline void Connectivity::init(
    const std::string& fname,
    const std::string& dataset,
    ElementType type,
    const std::vector<size_t>& shape)
{
    m_type = type;
    m_fname = fname;
    m_dataset = dataset;
    m_shape = shape;

    GOOSEFEM_ASSERT(m_shape.size() == 2);

#ifdef GOOSEFEM_ENABLE_ASSERT
    if (m_type == ElementType::Triangle) {
        GOOSEFEM_ASSERT(m_shape[1] == 3);
    }
    else if (m_type == ElementType::Quadrilateral) {
        GOOSEFEM_ASSERT(m_shape[1] == 4);
    }
    else if (m_type == ElementType::Hexahedron) {
        GOOSEFEM_ASSERT(m_shape[1] == 8);
    }
#endif
}

inline size_t Connectivity::nelem() const
{
    return m_shape[0];
}

inline size_t Connectivity::nne() const
{
    return m_shape[1];
}

inline std::vector<size_t> Connectivity::shape() const
{
    return m_shape;
}

inline std::string Connectivity::fname() const
{
    return m_fname;
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void Connectivity::checkShape()
{
    if (m_verified)
        return;

    H5Easy::File data(m_fname, H5Easy::File::ReadOnly);

    GOOSEFEM_CHECK(m_shape == H5Easy::getShape(data, m_dataset));

    m_verified = true;
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void Connectivity::readShape(const H5Easy::File& data)
{
    if (m_fname.size() == 0)
        m_fname = data.getName();

    m_shape = H5Easy::getShape(data, m_dataset);

    GOOSEFEM_ASSERT(m_shape.size() == 2);

#ifdef GOOSEFEM_ENABLE_ASSERT
    if (m_type == ElementType::Triangle) {
        GOOSEFEM_ASSERT(m_shape[1] == 3);
    }
    else if (m_type == ElementType::Quadrilateral) {
        GOOSEFEM_ASSERT(m_shape[1] == 4);
    }
    else if (m_type == ElementType::Hexahedron) {
        GOOSEFEM_ASSERT(m_shape[1] == 8);
    }
#endif

    m_verified = true;
}
#endif

inline std::vector<std::string> Connectivity::xdmf(size_t n_indent) const
{
    std::vector<std::string> ret;

    ret.push_back(
        "<Topology NumberOfElements=\"" + std::to_string(m_shape[0]) + "\" TopologyType=\"" +
        to_string(m_type) + "\">");

    ret.push_back(
        indent(n_indent) + "<DataItem Dimensions=\"" + std::to_string(m_shape[0]) + " " +
        std::to_string(m_shape[1]) + "\" Format=\"HDF\">" + m_fname + ":" + m_dataset +
        "</DataItem>");

    ret.push_back("</Topology>");

    return ret;
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Coordinates::Coordinates(const H5Easy::File& data, const std::string& dataset)
    : m_dataset(dataset)
{
    this->readShape(data);
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Coordinates::Coordinates(const std::string& fname, const std::string& dataset)
    : m_fname(fname), m_dataset(dataset)
{
    H5Easy::File data(m_fname, H5Easy::File::ReadOnly);

    this->readShape(data);
}
#endif

inline Coordinates::Coordinates(
    const std::string& fname, const std::string& dataset, const std::vector<size_t>& shape)
    : m_fname(fname), m_dataset(dataset), m_shape(shape)
{
    GOOSEFEM_ASSERT(m_shape.size() == 2);
}

inline size_t Coordinates::nnode() const
{
    return m_shape[0];
}

inline size_t Coordinates::ndim() const
{
    return m_shape[1];
}

inline std::vector<size_t> Coordinates::shape() const
{
    return m_shape;
}

inline std::string Coordinates::fname() const
{
    return m_fname;
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void Coordinates::checkShape()
{
    if (m_verified)
        return;

    H5Easy::File data(m_fname, H5Easy::File::ReadOnly);

    GOOSEFEM_CHECK(m_shape == H5Easy::getShape(data, m_dataset));

    m_verified = true;
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void Coordinates::readShape(const H5Easy::File& data)
{
    if (m_fname.size() == 0)
        m_fname = data.getName();

    m_shape = H5Easy::getShape(data, m_dataset);

    GOOSEFEM_ASSERT(m_shape.size() == 2);

    m_verified = true;
}
#endif

inline std::vector<std::string> Coordinates::xdmf(size_t n_indent) const
{
    std::vector<std::string> ret;

    if (m_shape[1] == 1) {
        ret.push_back("<Geometry GeometryType=\"X\">");
    }
    else if (m_shape[1] == 2) {
        ret.push_back("<Geometry GeometryType=\"XY\">");
    }
    else if (m_shape[1] == 3) {
        ret.push_back("<Geometry GeometryType=\"XYZ\">");
    }

    ret.push_back(
        indent(n_indent) + "<DataItem Dimensions=\"" + std::to_string(m_shape[0]) + " " +
        std::to_string(m_shape[1]) + "\" Format=\"HDF\">" + m_fname + ":" + m_dataset +
        "</DataItem>");

    ret.push_back("</Geometry>)");

    return ret;
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Attribute::Attribute(
    const H5Easy::File& data,
    const std::string& dataset,
    const std::string& name,
    AttributeType type)
    : m_type(type), m_dataset(dataset), m_name(name)
{
    this->readShape(data);
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline Attribute::Attribute(
    const std::string& fname,
    const std::string& dataset,
    const std::string& name,
    AttributeType type)
    : m_type(type), m_fname(fname), m_dataset(dataset), m_name(name)
{
    H5Easy::File data(m_fname, H5Easy::File::ReadOnly);

    this->readShape(data);
}
#endif

inline Attribute::Attribute(
    const std::string& fname,
    const std::string& dataset,
    const std::string& name,
    AttributeType type,
    const std::vector<size_t>& shape)
    : m_type(type), m_fname(fname), m_dataset(dataset), m_name(name), m_shape(shape)
{
    GOOSEFEM_ASSERT(m_shape.size() > 0);
}

inline std::vector<size_t> Attribute::shape() const
{
    return m_shape;
}

inline std::string Attribute::fname() const
{
    return m_fname;
}

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void Attribute::checkShape()
{
    if (m_verified)
        return;

    H5Easy::File data(m_fname, H5Easy::File::ReadOnly);

    GOOSEFEM_CHECK(m_shape == H5Easy::getShape(data, m_dataset));

    m_verified = true;
}
#endif

#ifndef GOOSEFEM_NO_HIGHFIVE
inline void Attribute::readShape(const H5Easy::File& data)
{
    if (m_fname.size() == 0) {
        m_fname = data.getName();
    }

    m_shape = H5Easy::getShape(data, m_dataset);

    GOOSEFEM_ASSERT(m_shape.size() > 0);

    m_verified = true;
}
#endif

inline std::vector<std::string> Attribute::xdmf(size_t n_indent) const
{
    GOOSEFEM_ASSERT(m_shape.size() > 0);
    GOOSEFEM_ASSERT(m_shape.size() < 3);

    std::vector<std::string> ret;

    if (m_shape.size() == 1) {
        ret.push_back(
            "<Attribute AttributeType=\"Scalar\" Center=\"" + to_string(m_type) + "\" Name=\"" +
            m_name + "\">");

        ret.push_back(
            indent(n_indent) + "<DataItem Dimensions=\"" + std::to_string(m_shape[0]) +
            "\" Format=\"HDF\">" + m_fname + ":" + m_dataset + "</DataItem>");
    }
    else if (m_shape.size() == 2) {
        ret.push_back(
            "<Attribute AttributeType=\"Vector\" Center=\"" + to_string(m_type) + "\" Name=\"" +
            m_name + "\">");

        ret.push_back(
            indent(n_indent) + "<DataItem Dimensions=\"" + std::to_string(m_shape[0]) + " " +
            std::to_string(m_shape[1]) + "\" Format=\"HDF\">" + m_fname + ":" + m_dataset +
            "</DataItem>");
    }

    ret.push_back("</Attribute>)");

    return ret;
}

inline Mesh::Mesh(const Connectivity& conn, const Coordinates& coor) : m_conn(conn), m_coor(coor)
{
}

inline void Mesh::push_back(const Attribute& data)
{
    m_attr.push_back(data);
}

inline std::vector<std::string> Mesh::xdmf(size_t n_indent) const
{
    std::vector<std::string> ret;

    {
        std::vector<std::string> lines = m_conn.xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(line);
        }
    }

    {
        std::vector<std::string> lines = m_coor.xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(line);
        }
    }

    for (auto& i : m_attr) {
        std::vector<std::string> lines = i.xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(line);
        }
    }

    return ret;
}

inline void Mesh::write(const std::string& fname, size_t n_indent) const
{
    TimeSeries({Increment(m_conn, m_coor, m_attr)}).write(fname, n_indent);
}

inline Increment::Increment(const Connectivity& conn, const Coordinates& coor)
{
    m_conn.push_back(conn);
    m_coor.push_back(coor);
}

inline Increment::Increment(
    const Connectivity& conn, const Coordinates& coor, const std::vector<Attribute>& attr)
    : m_attr(attr)
{
    m_conn.push_back(conn);
    m_coor.push_back(coor);
}

inline void Increment::push_back(const Connectivity& data)
{
    m_conn.push_back(data);
}

inline void Increment::push_back(const Coordinates& data)
{
    m_coor.push_back(data);
}

inline void Increment::push_back(const Attribute& data)
{
    m_attr.push_back(data);
}

inline std::vector<std::string> Increment::xdmf(size_t n_indent) const
{
    std::vector<std::string> ret;

    for (auto& i : m_conn) {
        std::vector<std::string> lines = i.xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(line);
        }
    }

    for (auto& i : m_coor) {
        std::vector<std::string> lines = i.xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(line);
        }
    }

    for (auto& i : m_attr) {
        std::vector<std::string> lines = i.xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(line);
        }
    }

    return ret;
}

inline TimeSeries::TimeSeries(const std::vector<Increment>& data) : m_data(data)
{
}

inline void TimeSeries::push_back(const Increment& data)
{
    m_data.push_back(data);
}

inline std::vector<std::string> TimeSeries::xdmf(size_t n_indent) const
{
    std::vector<std::string> ret;

    for (size_t inc = 0; inc < m_data.size(); ++inc) {
        ret.push_back("<Grid Name=\"Increment " + std::to_string(inc) + "\">");

        ret.push_back(indent(n_indent) + "<Time Value=\"" + std::to_string(inc) + "\"/>");

        std::vector<std::string> lines = m_data[inc].xdmf(n_indent);

        for (auto& line : lines) {
            ret.push_back(indent(n_indent) + line);
        }

        ret.push_back("</Grid>)");
    }

    return ret;
}

inline void TimeSeries::write(const std::string& fname, size_t n_indent) const
{
    std::ofstream myfile;

    myfile.open(fname);

    myfile << "<Xdmf Version=\"2.0\">" << std::endl;
    myfile << indent(n_indent) + "<Domain>" << std::endl;
    myfile << indent(n_indent * 2) +
                  "<Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"TimeSeries\">"
           << std::endl;

    std::vector<std::string> lines = this->xdmf(n_indent);

    for (auto& line : lines) {
        myfile << indent(n_indent * 3) + line << std::endl;
    }

    myfile << indent(n_indent * 2) + "</Grid>" << std::endl;
    myfile << indent(n_indent) + "</Domain>" << std::endl;
    myfile << "</Xdmf>" << std::endl;

    myfile.close();
}

} // namespace HDF5
} // namespace ParaView
} // namespace GooseFEM

#endif
