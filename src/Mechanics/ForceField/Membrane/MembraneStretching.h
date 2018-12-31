#ifndef MEDYAN_MembraneStretching_h
#define MEDYAN_MembraneStretching_h

#include "MembraneInteractions.h"

template <class MembraneStretchingInteractionType>
class MembraneStretching: public MembraneInteractions {
private:
    MembraneStretchingInteractionType _FFType;

public:
    virtual double computeEnergy(bool stretched) override;
    virtual void computeForces();
    virtual void computeForcesAux();

    virtual const string getName() { return "Membrane Stretching"; }
    
};



#endif
