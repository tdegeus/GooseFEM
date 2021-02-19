/**
Version information.

\file Version.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VERSION_H
#define GOOSEFEM_VERSION_H

#include "config.h"

namespace GooseFEM {

/**
Return git branch and hash (at the time of configuring).

\return ``{branch, hash}``
*/
inline std::vector<std::string> git();

/**
Return version as string.

\return Version string, e.g. ``v0.5.2``.
*/
inline std::string version();

/**
Return versions of this library and of all of its dependencies.

\return List of strings with version information.
*/
inline std::vector<std::string> version_dependencies();

} // namespace GooseFEM

#include "Version.hpp"

#endif
