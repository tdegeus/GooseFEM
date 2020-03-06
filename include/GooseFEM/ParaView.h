/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_PARAVIEW_H
#define GOOSEFEM_PARAVIEW_H

// See: http://xdmf.org/index.php/XDMF_Model_and_Format

#include "config.h"

#include <fstream>

#ifndef GOOSEFEM_NO_HIGHFIVE
#include <highfive/H5Easy.hpp>
#endif

namespace GooseFEM {
namespace ParaView {
namespace HDF5 {

inline std::string join(const std::vector<std::string>& lines, const std::string& sep = "\n");
inline std::string indent(size_t n);
xt::xtensor<double,2> as3d(const xt::xtensor<double,2>& data);

enum class ElementType {
    Triangle,
    Quadrilateral,
    Hexahedron };

enum class AttributeType {
    Cell,
    Node };

inline ElementType convert(GooseFEM::Mesh::ElementType type);
inline std::string to_string(ElementType type);
inline std::string to_string(AttributeType type);

class Connectivity {
public:
    Connectivity() = default;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Connectivity(
        const H5Easy::File& data, const std::string& dataset, GooseFEM::Mesh::ElementType type);
    #endif

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Connectivity(
        const std::string& fname, const std::string& dataset, GooseFEM::Mesh::ElementType type);
    #endif

    Connectivity(
        const std::string& fname,
        const std::string& dataset,
        GooseFEM::Mesh::ElementType type,
        const std::vector<size_t>& shape);

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Connectivity(const H5Easy::File& data, const std::string& dataset, ElementType type);
    #endif

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Connectivity(const std::string& fname, const std::string& dataset, ElementType type);
    #endif

    Connectivity(
        const std::string& fname,
        const std::string& dataset,
        ElementType type,
        const std::vector<size_t>& shape);

    size_t nelem() const;
    size_t nne() const;
    std::vector<size_t> shape() const;
    std::string fname() const;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    void checkShape();
    #endif

    std::vector<std::string> xdmf(size_t indent = 4) const;

private:
    #ifndef GOOSEFEM_NO_HIGHFIVE
    void readShape(const H5Easy::File& data);
    #endif

    #ifndef GOOSEFEM_NO_HIGHFIVE
    void init(const H5Easy::File& data, const std::string& dataset, ElementType type);
    #endif

    #ifndef GOOSEFEM_NO_HIGHFIVE
    void init(const std::string& fname, const std::string& dataset, ElementType type);
    #endif

    void init(
        const std::string& fname,
        const std::string& dataset,
        ElementType type,
        const std::vector<size_t>& shape);

    ElementType m_type;
    std::string m_fname;
    std::string m_dataset;
    std::vector<size_t> m_shape;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    bool m_verified = false; // if true: shape read from file, not need to check
    #endif
};

class Coordinates {
public:
    Coordinates() = default;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Coordinates(const H5Easy::File& data, const std::string& dataset);
    #endif

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Coordinates(
        const std::string& fname,
        const std::string& dataset);
    #endif

    Coordinates(
        const std::string& fname,
        const std::string& dataset,
        const std::vector<size_t>& shape);

    size_t nnode() const;
    size_t ndim() const;
    std::vector<size_t> shape() const;
    std::string fname() const;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    void checkShape();
    #endif

    std::vector<std::string> xdmf(size_t indent = 4) const;

private:
    #ifndef GOOSEFEM_NO_HIGHFIVE
    void readShape(const H5Easy::File& data);
    #endif

    std::string m_fname;
    std::string m_dataset;
    std::vector<size_t> m_shape;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    bool m_verified = false; // if true: shape read from file, not need to check
    #endif
};

class Attribute {
public:
    Attribute() = default;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Attribute(
        const H5Easy::File& data,
        const std::string& dataset,
        const std::string& name,
        AttributeType type);
    #endif

    #ifndef GOOSEFEM_NO_HIGHFIVE
    Attribute(
        const std::string& fname,
        const std::string& dataset,
        const std::string& name,
        AttributeType type);
    #endif

    Attribute(
        const std::string& fname,
        const std::string& dataset,
        const std::string& name,
        AttributeType type,
        const std::vector<size_t>& shape);

    std::vector<size_t> shape() const;
    std::string fname() const;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    void checkShape();
    #endif

    std::vector<std::string> xdmf(size_t indent = 4) const;

private:
    #ifndef GOOSEFEM_NO_HIGHFIVE
    void readShape(const H5Easy::File& data);
    #endif

    AttributeType m_type;
    std::string m_fname;
    std::string m_dataset;
    std::string m_name;
    std::vector<size_t> m_shape;

    #ifndef GOOSEFEM_NO_HIGHFIVE
    bool m_verified = false; // if true: shape read from file, not need to check
    #endif
};

class Mesh {
public:
    Mesh() = default;
    Mesh(const Connectivity& conn, const Coordinates& coor);

    void push_back(const Attribute& data);

    std::vector<std::string> xdmf(size_t indent = 4) const;

    void write(const std::string& fname, size_t n_indent = 4) const;

private:
    Connectivity m_conn;
    Coordinates m_coor;
    std::vector<Attribute> m_attr;
};

class Increment {
public:
    Increment() = default;

    Increment(const Connectivity& conn, const Coordinates& coor);

    Increment(
        const Connectivity& conn,
        const Coordinates& coor,
        const std::vector<Attribute>& attr);

    void push_back(const Connectivity& data);
    void push_back(const Coordinates& data);
    void push_back(const Attribute& data);

    std::vector<std::string> xdmf(size_t indent = 4) const;

private:
    std::vector<Connectivity> m_conn;
    std::vector<Coordinates> m_coor;
    std::vector<Attribute> m_attr;
};

class TimeSeries {
public:
    TimeSeries() = default;
    TimeSeries(const std::vector<Increment>& data);

    void push_back(const Increment& data);

    std::vector<std::string> xdmf(size_t indent = 4) const;

    void write(const std::string& fname, size_t n_indent = 4) const;

private:
    std::vector<Increment> m_data;
};

} // namespace HDF5
} // namespace ParaView
} // namespace GooseFEM

#include "ParaView.hpp"

#endif
