/**
\file version.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VERSION_HPP
#define GOOSEFEM_VERSION_HPP

#include "version.h"

namespace GooseFEM {

namespace detail {

    inline std::string unquote(const std::string& arg)
    {
        std::string ret = arg;
        ret.erase(std::remove(ret.begin(), ret.end(), '\"'), ret.end());
        return ret;
    }

}

inline std::string version()
{
    return detail::unquote(std::string(QUOTE(GOOSEFEM_VERSION)));
}

inline std::vector<std::string> version_dependencies()
{
    std::vector<std::string> ret;

    ret.push_back("goosefem=" + version());

    ret.push_back("xtensor=" +
        detail::unquote(std::string(QUOTE(XTENSOR_VERSION_MAJOR))) + "." +
        detail::unquote(std::string(QUOTE(XTENSOR_VERSION_MINOR))) + "." +
        detail::unquote(std::string(QUOTE(XTENSOR_VERSION_PATCH))));

    #if defined(GOOSEFEM_EIGEN) || defined(EIGEN_WORLD_VERSION)

        ret.push_back("eigen=" +
            detail::unquote(std::string(QUOTE(EIGEN_WORLD_VERSION))) + "." +
            detail::unquote(std::string(QUOTE(EIGEN_MAJOR_VERSION))) + "." +
            detail::unquote(std::string(QUOTE(EIGEN_MINOR_VERSION))));

    #endif

    return ret;
}

} // namespace GooseFEM

#endif
