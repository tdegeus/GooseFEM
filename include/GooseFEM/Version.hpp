/**
\file Version.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VERSION_HPP
#define GOOSEFEM_VERSION_HPP

#include "Version.h"

namespace GooseFEM {

inline std::vector<std::string> git()
{
    return std::vector<std::string>{std::string(GOOSEFEM_GIT_BRANCH),
                                    std::string(GOOSEFEM_GIT_HASH)};
}

inline std::string version()
{
    return std::to_string(GOOSEFEM_VERSION_MAJOR) + "." +
           std::to_string(GOOSEFEM_VERSION_MINOR) + "." +
           std::to_string(GOOSEFEM_VERSION_PATCH);
}

inline std::vector<std::string> version_dependencies()
{
    std::vector<std::string> ret;

    ret.push_back("goosefem=" +
        std::to_string(GOOSEFEM_VERSION_MAJOR) + "." +
        std::to_string(GOOSEFEM_VERSION_MINOR) + "." +
        std::to_string(GOOSEFEM_VERSION_PATCH));

    ret.push_back("xtensor=" +
        std::to_string(XTENSOR_VERSION_MAJOR) + "." +
        std::to_string(XTENSOR_VERSION_MINOR) + "." +
        std::to_string(XTENSOR_VERSION_PATCH));

    #if defined(GOOSEFEM_EIGEN) || defined(EIGEN_WORLD_VERSION)

        ret.push_back("eigen=" +
            std::to_string(EIGEN_WORLD_VERSION) + "." +
            std::to_string(EIGEN_MAJOR_VERSION) + "." +
            std::to_string(EIGEN_MINOR_VERSION));

    #endif

    return ret;
}


} // namespace GooseFEM

#endif
