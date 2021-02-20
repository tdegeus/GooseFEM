/**
Version information.

\file version.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VERSION_H
#define GOOSEFEM_VERSION_H

#include "config.h"

/**
Current version.

Either:

-   Configure using CMake at install time. Internally uses::

        python -c "from setuptools_scm import get_version; print(get_version())"

-   Define externally using::

        -DGOOSEFEM_VERSION="`python -c "from setuptools_scm import get_version; print(get_version())"`"

    From the root of this project. This is what ``setup.py`` does.

Note that both ``CMakeLists.txt`` and ``setup.py`` will construct the version string using
*setuptools_scm* **unless** an environment ``PKG_VERSION`` is defined.
If ``PKG_VERSION`` is defined the version string will be read from that variable.
*/
#ifndef GOOSEFEM_VERSION
#define GOOSEFEM_VERSION "@GOOSEFEM_VERSION@"
#endif

namespace GooseFEM {

/**
Return version string, e.g. "0.8.0"

\return std::string
*/
inline std::string version();

/**
Return versions of this library and of all of its dependencies.
The output is a list of strings, e.g.::

    "goosefem=0.7.0",
    "xtensor=0.20.1"
    ...

\return List of strings.
*/
inline std::vector<std::string> version_dependencies();

} // namespace GooseFEM

#include "version.hpp"

#endif
