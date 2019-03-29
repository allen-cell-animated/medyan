#ifndef MEDYAN_MembraneStretching_h
#define MEDYAN_MembraneStretching_h

#include "MembraneInteractions.h"

template <class MembraneStretchingInteractionType>
class MembraneStretching: public MembraneInteractions {
private:
    MembraneStretchingInteractionType _FFType;

public:
    virtual double computeEnergy(const double* coord, bool stretched) override;
    virtual void computeForces(const double* coord, double* force) override;
    virtual void computeForcesAux();

    virtual const string getName() { return "Membrane Stretching"; }
    
};



#endif
