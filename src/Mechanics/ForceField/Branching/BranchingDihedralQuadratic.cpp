#include "Mechanics/ForceField/Branching/BranchingDihedralQuadratic.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint> // uint_fast8_t
#include <limits>

#include "Mechanics/ForceField/Branching/BranchingDihedral.h"
#include "Structure/BranchingPoint.h"
#include "Util/Io/Log.hpp"
#include "Util/Math/Vec.hpp"

using namespace mathfunc;
using V3 = Vec< 3, floatingpoint >;

constexpr floatingpoint cosDihTol = 0.01;
constexpr floatingpoint angSinMin = 0.001;
constexpr floatingpoint dihSinMin = 0.00001;

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
        const auto p = pos[i];
        const auto mp = (1 - p) * coord1 + p * coord2;

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
        auto ct = dot(n1, n2) / (n1mag * n2mag);

        if(ct > 1.0 + cosDihTol || ct < -1.0 - cosDihTol) {
            LOG(WARNING) << "Dihedral incorrect cos theta: " << ct;
            LOG(INFO) << "Interaction information:\n"
                << "Bead coords: " << coord1 << ' ' << coord2 << ' ' << coord3 << ' ' << coord4 << '\n'
                << "Position: " << p << ", brancher coord: " << mp;
        }

        if(ct < -1.0) ct = -1.0;
        if(ct >  1.0) ct =  1.0;

        const auto theta = std::acos(ct);

        const auto U_i = kdih[i] * theta * theta;

        if(!std::isfinite(U_i) || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return std::numeric_limits<floatingpoint>::infinity();
        }

        U += U_i;
    } // End loop interactions

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

    for(size_t i = 0; i < nint; ++i) {
        auto coord1 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i    ]);
        auto coord2 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i + 1]);
        auto coord3 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i + 2]);
        auto coord4 = makeRefVec< 3, floatingpoint >(coord + 3 * beadSet[bpi * i + 3]);

        auto f1 = makeRefVec< 3, floatingpoint >(f + 3 * beadSet[bpi * i    ]);
        auto f2 = makeRefVec< 3, floatingpoint >(f + 3 * beadSet[bpi * i + 1]);
        auto f3 = makeRefVec< 3, floatingpoint >(f + 3 * beadSet[bpi * i + 2]);
        auto f4 = makeRefVec< 3, floatingpoint >(f + 3 * beadSet[bpi * i + 3]);

        // Brancher coordinate on the mother filament
        const auto p = pos[i];
        const auto mp = (1 - p) * coord1 + p * coord2;

        // Bonds
        const auto b1 = coord2 - mp;
        const auto b2 = coord3 - mp;
        const auto b3 = coord4 - coord3;

        //---------------------------------------------------------------------
        // E = E1(b1, b2, b3)
        //   = E2(c2, mp, c3, c4)
        //   = E3(c1, c2, c3, c4)
        //
        // If we have dE1, and let E1i = dE1 / dbi, then
        //   dE2 / dc2 = E11
        //   dE2 / dmp = - E11 - E12
        //   dE2 / dc3 = E12 - E13
        //   dE2 / dc4 = E13
        //
        //   dE3 / dc1 = - (1-p)E11 - (1-p)E12
        //   dE3 / dc2 = (1-p)E11 - p E12
        //   dE3 / dc3 = E12 - E13
        //   dE3 / dc4 = E13
        //
        //---------------------------------------------------------------------
        //                 (b1 x b2) . (b3 x b2)
        // cos(theta) = ----------------------------
        //                 |b1 x b2|   |b3 x b2|
        //
        //                 (b1 . b3) (b2 . b2) - (b1 . b2) (b3 . b2)
        //            = -----------------------------------------------
        //                 |b1| |b3| |b2|^2  sin<b1, b2> sin<b3, b2>
        //
        //                 cos<b1, b3> - cos<b1, b2> cos<b3, b2>
        //            = -------------------------------------------
        //                       sin<b1, b2>  sin<b3, b2>
        //
        //---------------------------------------------------------------------
        // Useful formula
        //
        // d cos<x1, x2> / dx1
        //
        //           x2         (x1 . x2)  x1
        //    = ----------- - -----------------
        //       |x1| |x2|       |x1|^3 |x2|
        //
        //                  [    x2          x1    ]
        //    = cos<x1, x2> [ --------- - -------- ]
        //                  [  x1 . x2     |x1|^2  ]
        //
        // d sin<x1, x2> / dx1
        //
        //         cos<x1, x2>
        //    = - ------------- d cos<x1, x2> / dx1
        //         sin<x1, x2>
        //
        //---------------------------------------------------------------------
        const auto b1mag2 = magnitude2(b1);
        const auto b2mag2 = magnitude2(b2);
        const auto b3mag2 = magnitude2(b3);
        const auto b1mag  = std::sqrt(b1mag2);
        const auto b2mag  = std::sqrt(b2mag2);
        const auto b3mag  = std::sqrt(b3mag2);

        const auto mag13inv = (floatingpoint)1.0 / (b1mag * b3mag);
        const auto mag12inv = (floatingpoint)1.0 / (b1mag * b2mag);
        const auto mag32inv = (floatingpoint)1.0 / (b3mag * b2mag);

        const auto dot13 = dot(b1, b3);
        const auto dot12 = dot(b1, b2);
        const auto dot32 = dot(b3, b2);

        // Bond angles
        const auto cos12     = dot12 * mag12inv;
        const auto sin12_2   = std::max< floatingpoint >(1 - cos12 * cos12, 0.0);
        const auto sin12     = std::max< floatingpoint >(std::sqrt(sin12_2), angSinMin);
        const auto sin12inv  = (floatingpoint)1.0 / sin12;
        const auto sin12inv2 = sin12inv * sin12inv;

        const auto cos32     = dot32 * mag32inv;
        const auto sin32_2   = std::max< floatingpoint >(1 - cos32 * cos32, 0.0);
        const auto sin32     = std::max< floatingpoint >(std::sqrt(sin32_2), angSinMin);
        const auto sin32inv  = (floatingpoint)1.0 / sin32;
        const auto sin32inv2 = sin32inv * sin32inv;

        const auto cos13 = dot13 * mag13inv;

        // ct -- cos(theta)
        const auto ctFac1 = cos13 - cos12 * cos32; // ct numerator
        const auto ctFac2 = sin12inv * sin32inv;   // ct inv denominator
        const auto ct     = ctFac1 * ctFac2;

        if(ct > 1.0 + cosDihTol || ct < -1.0 - cosDihTol) {
            LOG(WARNING) << "Dihedral incorrect cos theta: " << ct;
            LOG(INFO) << "Interaction information:\n"
                << "Bead coords: " << coord1 << ' ' << coord2 << ' ' << coord3 << ' ' << coord4 << '\n'
                << "Position: " << p << ", brancher coord: " << mp;
        }

        if(ct < -1.0) ct = -1.0;
        if(ct >  1.0) ct =  1.0;

        // derivatives on ct
        //---------------------------------------------------------------------
        // d(ct_nu) / db1
        //
        //           b3                      b1
        //    = ----------- - cos<b1, b3> --------
        //       |b1| |b3|                 |b1|^2
        //
        //                        b2                                   b1
        //      - cos<b3, b2> ----------- + cos<b3, b2> cos<b1, b2> --------
        //                     |b1| |b2|                             |b1|^2
        //
        // d(ct_nu) / db2
        //
        //                         b3                                  b2
        //    = - cos<b1, b2> ----------- + cos<b1, b2> cos<b3, b2> --------
        //                     |b2| |b3|                             |b2|^2
        //
        //                         b1                                  b2
        //      - cos<b3, b2> ----------- + cos<b1, b2> cos<b3, b2> --------
        //                     |b2| |b1|                             |b2|^2
        //
        // d(ln ct_de) / db1
        //
        //         cos^2 <b1, b2>  [     b2         b1    ]
        //    = - ---------------- [ --------- - -------- ]
        //         sin^2 <b1, b2>  [  b1 . b2     |b1|^2  ]
        //
        // d(ln ct_de) / db2
        //
        //         cos^2 <b3, b2>  [     b3         b2    ]
        //    = - ---------------- [ --------- - -------- ]
        //         sin^2 <b3, b2>  [  b3 . b2     |b2|^2  ]
        //
        //         cos^2 <b1, b2>  [     b1         b2    ]
        //      - ---------------- [ --------- - -------- ]
        //         sin^2 <b1, b2>  [  b1 . b2     |b2|^2  ]
        //
        //---------------------------------------------------------------------
        //          d(ct_nu)
        // d(ct) = ---------- - ct * d(ln ct_de)
        //            ct_de
        //
        // Let E1i = E1i1 * b1 + E1i2 * b2 + E1i3 * b3, then
        //
        //             ct       cos^2 <b1, b2>       ct                  ct
        // E111 = - -------- - ---------------- * -------- = - -----------------------
        //           |b1|^2     sin^2 <b1, b2>     |b1|^2       sin^2 <b1, b2> |b1|^2
        //
        //             cos<b3, b2>       cos^2 <b1, b2>   ct
        // E112 = - ----------------- + -------------------------
        //           ct_de |b1| |b2|     sin^2 <b1, b2>  b1 . b2
        //
        //             1      [    cos<b3, b2>     cos<b1, b2> ct  ]
        //      = ----------- [ - ------------- + ---------------- ]
        //         |b1| |b2|  [       ct_de        sin^2 <b1, b2>  ]
        //
        // E113 = 1 / (|b1| |b3| ct_de)
        //
        //             cos<b3, b2>       cos^2 <b1, b2>   ct
        // E121 = - ----------------- + ------------------------- = E112
        //           ct_de |b1| |b2|     sin^2 <b1, b2>  b1 . b2
        //
        //         2 cos<b1, b2> cos<b3, b2>      cos^2 <b1, b2>  ct        cos^2 <b3, b2>  ct
        // E122 = --------------------------- - ----------------------- - -----------------------
        //              ct_de   |b2|^2           sin^2 <b1, b2> |b2|^2     sin^2 <b3, b2> |b2|^2
        //
        //            1    [  2 cos<b1, b3>       [      cos^2 <b1, b2>     cos^2 <b3, b2>  ] ]
        //      = -------- [ --------------- - ct [ 2 + ---------------- + ---------------- ] ]
        //         |b2|^2  [      ct_de           [      sin^2 <b1, b2>     sin^2 <b3, b2>  ] ]
        //
        //            1    [  2 cos<b1, b3>       [         1                  1        ] ]
        //      = -------- [ --------------- - ct [ ---------------- + ---------------- ] ]
        //         |b2|^2  [      ct_de           [  sin^2 <b1, b2>     sin^2 <b3, b2>  ] ]
        //
        //             cos<b1, b2>       cos^2 <b3, b2>   ct
        // E123 = - ----------------- + ------------------------- = E132
        //           ct_de |b3| |b2|     sin^2 <b3, b2>  b3 . b2
        //
        //---------------------------------------------------------------------

    // TODO
    } // End loop interactions

} // void force(...)
