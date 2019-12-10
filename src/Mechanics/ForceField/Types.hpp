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

#endif