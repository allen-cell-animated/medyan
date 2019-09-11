#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_Hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_Hpp

#include "Mechanics/ForceField/Membrane/MembraneInteractions.hpp"

enum class MembraneStretchingAccumulationType {
    ByTriangle, ByVertex
};

template< typename Impl, MembraneStretchingAccumulationType >
class MembraneStretching: public MembraneInteractions {
private:
    Impl _impl;

public:
    virtual floatingpoint computeEnergy(const floatingpoint* coord, bool stretched) override;
    virtual void computeForces(const floatingpoint* coord, floatingpoint* force) override;

    virtual string getName() const override { return "Membrane Stretching"; }
    
};

// Specialization declarations
// template< typename Impl > double MembraneStretching< Impl, MembraneStretchingAccumulationType::ByTriangle >::computeEnergy(const double* coord, bool stretched);
// template< typename Impl > void MembraneStretching< Impl, MembraneStretchingAccumulationType::ByTriangle >::computeForces(const double* coord, double* force);
// template< typename Impl > double MembraneStretching< Impl, MembraneStretchingAccumulationType::ByVertex >::computeEnergy(const double* coord, bool stretched);
// template< typename Impl > void MembraneStretching< Impl, MembraneStretchingAccumulationType::ByVertex >::computeForces(const double* coord, double* force);


#endif
