/**
\file version.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VERSION_HPP
#define GOOSEFEM_VERSION_HPP

#include "version.h"

namespace GooseFEM {

inline std::string version()
{
    return std::string(QUOTE(GOOSEFEM_VERSION));
}

inline std::vector<std::string> version_dependencies()
{
    std::vector<std::string> ret;

    ret.push_back("goosefem=" + version());

    ret.push_back("xtensor=" +
        std::to_string(QUOTE(XTENSOR_VERSION_MAJOR)) + "." +
        std::to_string(QUOTE(XTENSOR_VERSION_MINOR)) + "." +
        std::to_string(QUOTE(XTENSOR_VERSION_PATCH)));

    #if defined(GOOSEFEM_EIGEN) || defined(EIGEN_WORLD_VERSION)

        ret.push_back("eigen=" +
            std::to_string(QUOTE(EIGEN_WORLD_VERSION)) + "." +
            std::to_string(QUOTE(EIGEN_MAJOR_VERSION)) + "." +
            std::to_string(QUOTE(EIGEN_MINOR_VERSION)));

    #endif

    return ret;
}

} // namespace GooseFEM

#endif
