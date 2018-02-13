#ifndef MEDYAN_VolumeConservationMembrane_h
#define MEDYAN_VolumeConservationMembrane_h

#include "VolumeConservationInteractions.h"

template<class VolumeConservationMembraneInteractionType>
class VolumeConservationMembrane: public VolumeConservationInteractions {
private:
    VolumeConservationMembraneInteractionType _FFType;

public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();

    virtual const string getName()const { return "Membrane Volume Conservation"; }
};


#endif
