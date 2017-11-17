#ifndef MEDYAN_MembraneBending_h
#define MEDYAN_MembraneBending_h

#include "common.h"

#include "MembraneInteractions.h"

//FORWARD DECLARATIONS
class Membrane;

/// Represents a Filament bending interaction
template <class MembraneBendingInteractionType>
class MembraneBending : public MembraneInteractions {
    
private:
    MembraneBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Membrane Bending";}
};

#endif
