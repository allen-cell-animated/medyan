#ifndef MEDYAN_VolumeConservationMembrane_h
#define MEDYAN_VolumeConservationMembrane_h

#include "Mechanics/ForceField/VolumeConservation/Interactions.hpp"

template<class InteractionType>
class VolumeConservationMembrane: public VolumeConservationInteractions {
private:
    InteractionType _FFType;

public:
    virtual floatingpoint computeEnergy(const floatingpoint* coord, bool stretched) override;
    virtual void computeForces(const floatingpoint* coord, floatingpoint* force) override;

    virtual std::string getName() const override { return "Membrane Volume Conservation"; }
};


#endif
