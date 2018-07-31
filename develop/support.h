
#ifndef SUPPORT_H
#define SUPPORT_H

#include <catch2/catch.hpp>

#include <cppmat/cppmat.h>

#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>

#include "../include/xGooseFEM/GooseFEM.h"

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );


template <class E, class W>
inline auto Average(E&& e, W&& weights)
{
    return xt::sum(e*weights, {0,1}) / xt::sum(weights, {0,1});
}


#endif
