
#ifndef SUPPORT_H
#define SUPPORT_H

#include <catch/catch.hpp>
#include <cppmat/cppmat.h>
#include <Eigen/Eigen>

#include "../src/GooseFEM/GooseFEM.h"

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#endif
