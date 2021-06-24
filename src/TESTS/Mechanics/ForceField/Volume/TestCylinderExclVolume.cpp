#include <numeric> // iota
#include <type_traits> // is_same
#include <vector>

#include "catch2/catch.hpp"

#include "Mechanics/ForceField/Volume/CylinderExclVolRepulsion.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

using namespace test_ff_common;

TEST_CASE("Force field: Cylinder excl volume", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace medyan;
    using VF = std::vector< floatingpoint >;

    // Prepare data
    //---------------------------------
    struct TestInput {
        std::string      name;
        std::vector<int> beadSet;
        VF               k;
        VF               coord;
    };


    // Multiple test sets
    std::vector< TestInput > testInputs;

    // Input: standard symmetric case
    testInputs.push_back({
        "standard-symmetric",
        // bead set
        { 0, 1, 2, 3 },
        // k
        { 100.0 },
        // coordinates
        {
            10, 0, 0,
            -10, 0, 0,

            0, 10, 5,
            0, -10, 5,
        },
    });

    // Input: general cases
    testInputs.push_back({
        "general",
        // bead set
        { 0, 1, 2, 3,  0, 1, 4, 5 },
        // k
        { 100.0, 2e3 },
        // coordinates
        {
            10, 0, 0,
            -10, 0, 0,

            0, 10, 5,
            3, -10, 5,

            3, 70, 30,
            2, -20, 0,
        },
    });

    // Prepare test parameters
    //---------------------------------
    const floatingpoint moveMag      = std::is_same< floatingpoint, float >::value ? 5e-3 : 1e-5;

    const floatingpoint diffDeRelEps = std::is_same< floatingpoint, float >::value ? 6e-2 : 5e-4;

    const size_t repsTot     = 10;
    const size_t repsPassReq = 9;


    // Run tests for each set of input
    //---------------------------------
    for(auto& ti : testInputs) {

        // Prepare functions
        //---------------------------------
        const auto calcEnergy = [&](const VF& c) {
            auto tempC = c;
            auto ret = CylinderExclVolRepulsion{}.energy(
                tempC.data(),
                ti.beadSet.data(),
                ti.k.data(),
                ti.k.size()
            );
            return ret;
        };
        const auto calcForce = [&](const VF& c, VF& f) {
            auto tempC = c;
            CylinderExclVolRepulsion{}.forces(
                tempC.data(),
                f.data(),
                ti.beadSet.data(),
                ti.k.data(),
                ti.k.size()
            );
        };


        // Run the test
        //---------------------------------
        size_t repsPass = 0;
        for(size_t rep = 1; rep <= repsTot; ++rep) {
            const auto res = testEnergyForceConsistency(ti.coord, calcEnergy, calcForce, moveMag, diffDeRelEps);
            if(res.passed)
                ++repsPass;
            else
                WARN(ti.name << " (Rep " << rep << '/' << repsTot << " fail) E: " << calcEnergy(ti.coord) << " Actual de: " << res.deActual << " Expected de: " << res.deExpected);
        }
        REQUIRE(repsPass >= repsPassReq);
    }
}
