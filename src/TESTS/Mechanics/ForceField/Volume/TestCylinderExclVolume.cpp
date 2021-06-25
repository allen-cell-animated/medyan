#include <numeric> // iota
#include <type_traits> // is_same
#include <vector>

#include "catch2/catch.hpp"

#include "Mechanics/ForceField/Volume/CylinderExclVolRepulsion.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

using namespace test_ff_common;

TEST_CASE("Force field: Cylinder excl volume scaling", "[ForceField]") {
    // The test case checks the scaling behavior of volume exclusion interaction

    using namespace std;

    const std::vector<floatingpoint> coord {
        // Segments of cylinder 1
        -60.0, 0.0, 0.0,
        0.0,   0.0, 0.0,
        60.0,  0.0, 0.0,

        // Segments of cylinder 2
        -40.0, -40.0, 10.0,
        -28.0, -24.0, 10.0,
        -10.0, 0.0,   10.0,
        20.0,  40.0,  10.0,
    };

    // Configuration 0
    const std::vector<int>           beadSet0  { 0, 2, 3, 6 };
    const std::vector<floatingpoint> k0        { 1e-2 };
    const std::vector<floatingpoint> eqLength0 { 100.0, 100.0 };

    // Configuration 1: breaking down each cylinder into several segments
    const std::vector<int>           beadSet1  { 0, 1, 3, 4,  0, 1, 4, 5,  0, 1, 5, 6,  1, 2, 3, 4,  1, 2, 4, 5,  1, 2, 5, 6, };
    const std::vector<floatingpoint> k1        { 1e-2,        1e-2,        1e-2,        1e-2,        1e-2,        1e-2,       };
    const std::vector<floatingpoint> eqLength1 { 50.0, 20.0,  50.0, 30.0,  50.0, 50.0,  50.0, 20.0,  50.0, 30.0,  50.0, 50.0, };

    auto tempCoord = coord;
    const auto energy0 = CylinderExclVolRepulsion{}.energy(tempCoord.data(), beadSet0.data(), k0.data(), eqLength0.data(), k0.size());
    const auto energy1 = CylinderExclVolRepulsion{}.energy(tempCoord.data(), beadSet1.data(), k1.data(), eqLength1.data(), k1.size());

    REQUIRE(energy0 == Approx(energy1));
}


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
        VF               eqLength;
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
        { 0.25 },
        // equilibrium length
        { 20.0, 20.0 },
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
        // equilibrium length
        { 20.0, 20.0,  20.0, 100.0 },
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
                ti.eqLength.data(),
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
                ti.eqLength.data(),
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
