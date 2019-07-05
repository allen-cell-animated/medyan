#ifndef MEDYAN_Mechanics_ForceField_Branching_BranchingDihedralQuadratic_Hpp
#define MEDYAN_Mechanics_ForceField_Branching_BranchingDihedralQuadratic_Hpp

#include "common.h" // floatingpoint

struct BranchingDihedralQuadratic {
    floatingpoint energy(
        const floatingpoint *coord, floatingpoint* f /* to be removed */, const int *beadSet,
        const floatingpoint *kdih, const floatingpoint *pos
    ) const;

    void forces(
        const floatingpoint *coord, floatingpoint *f, const int *beadSet,
        const floatingpoint *kdih, const floatingpoint *pos
    ) const;
};

#endif
