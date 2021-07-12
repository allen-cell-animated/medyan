#ifndef MEDYAN_Mechanics_ForceField_Types_Hpp
#define MEDYAN_Mechanics_ForceField_Types_Hpp

#include <vector>

#include "utility.h" // floatingpoint

struct ForceFieldTypes {
    enum class LoadForceEnd { Plus, Minus };
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

    // Total number of degrees of freedom. It is not necessarily the size of vectorized data, as the vectorized data may contain dependent variables.
    int ndof = 0;

};

#endif
