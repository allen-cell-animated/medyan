#include "Mechanics/ForceField/Branching/BranchingDihedralQuadratic.hpp"

#include <cmath> // acos, isfinite
#include <cstdint> // uint_fast8_t
#include <limits>

#include "Mechanics/ForceField/Branching/BranchingDihedral.h"
#include "Structure/BranchingPoint.h"
#include "Util/Math/Vec.hpp"

using namespace mathfunc;
using V3 = Vec< 3, floatingpoint >;

floatingpoint BranchingDihedralQuadratic::energy(
    const floatingpoint *coord, floatingpoint* f /* to be removed */, const int *beadSet,
    const floatingpoint *kdih, const floatingpoint *pos
) const {

    constexpr std::uint_fast8_t bpi = 4;
    static_assert(
        bpi == BranchingDihedral< BranchingDihedralQuadratic >::n,
        "Number of beads per interaction in branching dihedral quadratic does not match"
    );
    const auto nint = BranchingPoint::getBranchingPoints().size(); // should be passed as an argument

    floatingpoint U = 0.0;

    for(size_t i = 0; i < nint; ++i) {
        auto coord1 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i    ]);
        auto coord2 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i + 1]);
        auto coord3 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i + 2]);
        auto coord4 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i + 3]);

        // Brancher coordinate on the mother filament
        const auto mp = (1 - pos[i]) * coord1 + pos[i] * coord2;

        // Bonds
        const auto b1 = coord2 - mp;
        const auto b2 = coord3 - mp;
        const auto b3 = coord4 - coord3;

        // Triangle normals
        const auto n1    = cross(b1, b2);
        const auto n2    = cross(b3, b2);
        const auto n1mag = magnitude(n1);
        const auto n2mag = magnitude(n2);

        // cos(theta)
        auto c = dot(n1, n2) / (n1mag * n2mag);
        if(c < -1.0) c = -1.0;
        if(c >  1.0) c =  1.0;

        const auto theta = std::acos(c);

        const auto U_i = kdih[i] * theta * theta;

        if(!std::isfinite(U_i) || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return std::numeric_limits<floatingpoint>::infinity();
        }

        U += U_i;
    }

    return U;

} // floatingpoint energy(...)

void BranchingDihedralQuadratic::forces(
    const floatingpoint *coord, floatingpoint *f, const int *beadSet,
    const floatingpoint *kdih, const floatingpoint *pos
) const {

    constexpr std::uint_fast8_t bpi = 4;
    static_assert(
        bpi == BranchingDihedral< BranchingDihedralQuadratic >::n,
        "Number of beads per interaction in branching dihedral quadratic does not match"
    );
    auto nint = BranchingPoint::getBranchingPoints().size(); // should be passed as an argument

    // TODO

} // void force(...)
