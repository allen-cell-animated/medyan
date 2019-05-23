#ifndef MEDYAN_Mechancis_ForceField_Membrane_MembraneBending_Hpp
#define MEDYAN_Mechancis_ForceField_Membrane_MembraneBending_Hpp

#include "common.h"

#include "Mechanics/ForceField/Membrane/MembraneInteractions.hpp"

//FORWARD DECLARATIONS
class Membrane;

/// Represents a Filament bending interaction
template <class MembraneBendingInteractionType>
class MembraneBending : public MembraneInteractions {
    
private:
    MembraneBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(const double* coord, bool stretched) override;
    virtual void computeForces(const double* coord, double* force) override;
    
    virtual string getName() const override { return "Membrane Bending"; }
};

#endif
