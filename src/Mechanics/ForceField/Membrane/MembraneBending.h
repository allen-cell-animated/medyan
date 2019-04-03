#ifndef MEDYAN_MembraneBending_h
#define MEDYAN_MembraneBending_h

#include "common.h"

#include "Mechanics/ForceField/Membrane/MembraneInteractions.h"

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
    
    virtual const string getName() {return "Membrane Bending";}
};

#endif
