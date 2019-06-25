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
    virtual floatingpoint computeEnergy(const floatingpoint* coord, bool stretched) override;
    virtual void computeForces(const floatingpoint* coord, floatingpoint* force) override;
    
    virtual string getName() const override { return "Membrane Bending"; }
};

#endif
