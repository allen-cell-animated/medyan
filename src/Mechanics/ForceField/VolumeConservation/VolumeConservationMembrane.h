#ifndef MEDYAN_VolumeConservationMembrane_h
#define MEDYAN_VolumeConservationMembrane_h

#include "VolumeConservationInteractions.h"

template<class VolumeConservationMembraneInteractionType>
class VolumeConservationMembrane: public VolumeConservationInteractions {
private:
    VolumeConservationMembraneInteractionType _FFType;

public:
    virtual double computeEnergy(const double* coord, bool stretched) override;
    virtual void computeForces(const double* coord, double* force) override;
    virtual void computeForcesAux();

    virtual std::string getName()const { return "Membrane Volume Conservation"; }
};


#endif
