#ifndef MEDYAN_Mechanics_Minimizer_MinimizationTypes_Hpp
#define MEDYAN_Mechanics_Minimizer_MinimizationTypes_Hpp

#include "Mechanics/ForceField/Types.hpp"

struct MinimizationResult {

    EnergyReport energiesBefore;
    EnergyReport energiesAfter;
};

#endif
