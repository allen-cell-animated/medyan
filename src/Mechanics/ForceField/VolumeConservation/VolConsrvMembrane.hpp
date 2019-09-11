#ifndef MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvMembrane_hpp
#define MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvMembrane_hpp

#include "Mechanics/ForceField/VolumeConservation/VolConsrvInteractions.hpp"

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
