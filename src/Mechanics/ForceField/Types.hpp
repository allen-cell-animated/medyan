#ifndef MEDYAN_Mechanics_ForceField_Types_hpp
#define MEDYAN_Mechanics_ForceField_Types_hpp

// Notes
//   - This file is used to define various types about the force field, used by
//     various files. Do not include files here unless necessary, to avoid
//     circular includes.

#include <vector>

#include "utility.h" // floatingpoint

struct ForceFieldTypes {
    enum class LoadForceEnd { Plus, Minus };
    enum class GeometryCurvRequirement { curv, curv2 };
};

struct EnergyReport {
    struct EachEnergy {
        string name;
        floatingpoint energy;
    };

    floatingpoint total;
    std::vector< EachEnergy > individual;
};

// This information is used in the force field vectorization
struct FFCoordinateStartingIndex {
    std::size_t bead = 0;
    std::size_t bubble = 0;
    std::size_t vertex = 0;
    std::size_t mem2d = 0;
};

#endif
