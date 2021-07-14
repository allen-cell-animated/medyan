

#include "catch2/catch.hpp"

#include "Chemistry/ChemSimpleGillespieImpl.h"

TEST_CASE("ChemSimpleGillespieImpl tests", "[ChemSim]") {

    auto a= Species("A", 0, 1000000, REG);
    auto b= Species("B", 0, 1000000, REG);
    auto c= Species("C", 0, 1000000, REG);
    auto reaction= Reaction({&a,&b,&c}, 1.0);
    REQUIRE(1==2);
}