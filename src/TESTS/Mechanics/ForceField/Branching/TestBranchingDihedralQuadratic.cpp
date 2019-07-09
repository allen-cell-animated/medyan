#include <numeric> // iota
#include <vector>

#include "catch2/catch.hpp"

#include "Mechanics/ForceField/Branching/BranchingDihedralQuadratic.hpp"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

using namespace test_ff_common;

TEST_CASE("Force field: Branching Dihedral Quadratic", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace mathfunc;
    using V3 = Vec< 3, floatingpoint >;
    using VA3 = VecArray< 3, floatingpoint >;

    // Prepare data
    //---------------------------------
    const size_t bpi = 4;

    vector< floatingpoint > kdih;
    vector< floatingpoint > pos;
    VA3 coords;

    // small dihedral angle
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    coords.push_back(V3 { 1.5, 1.0, 0.1 });
    pos.push_back(0.6);
    kdih.push_back(1.0);

    // 90 deg dihedral angle (<90)
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    coords.push_back(V3 { 1.5, 0.1, 1.0 });
    pos.push_back(0.5);
    kdih.push_back(1.0);

    // 90 deg dihedral angle (>90)
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    coords.push_back(V3 { 1.5, -0.1, 1.0 });
    pos.push_back(0.5);
    kdih.push_back(1.0);

    // large dihedral angle
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    coords.push_back(V3 { 1.5, -1.0, 0.1 });
    pos.push_back(0.4);
    kdih.push_back(1.0);

    const auto nint = kdih.size();
    vector< unsigned > beadSet(coords.size());
    std::iota(beadSet.begin(), beadSet.end(), 0u);

    // Prepare functions
    //---------------------------------
    const auto calcEnergy = [&](const VA3& c) {
        return BranchingDihedralQuadratic {}.energy(
            c.data(), nint, beadSet.data(), kdih.data(), pos.data()
        );
    };
    const auto calcForce = [&](const VA3& c, VA3& f) {
        BranchingDihedralQuadratic {}.forces(
            c.data(), f.data(), nint, beadSet.data(), kdih.data(), pos.data()
        );
    };

    // Prepare test parameters
    //---------------------------------
    const floatingpoint moveMag = 1e-3;

    const floatingpoint diffDeRelEps = 1e-3;

    // Run the test
    //---------------------------------
    const auto res = testEnergyForceConsistency(coords, calcEnergy, calcForce, moveMag, diffDeRelEps);
    REQUIRE(res.passed);
}
