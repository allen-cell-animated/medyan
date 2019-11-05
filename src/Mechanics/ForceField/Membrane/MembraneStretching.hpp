#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp

#include "Mechanics/ForceField/ForceField.h"
#include "utility.h"

enum class MembraneStretchingAccumulationType {
    ByTriangle, ByVertex
};

template< typename Impl, MembraneStretchingAccumulationType >
class MembraneStretching: public ForceField {
private:
    Impl _impl;

    // Force buffer
    std::vector< floatingpoint > forceBuffer_;

public:
    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override;
    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override;

    virtual string getName() override { return "Membrane Stretching"; }

    const auto& getForceBuffer() const { return forceBuffer_; }

    // Useless overrides
    virtual void vectorize() override {}
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual void whoIsCulprit() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return std::vector<NeighborList*>(); }
};

// Specialization declarations
// template< typename Impl > double MembraneStretching< Impl, MembraneStretchingAccumulationType::ByTriangle >::computeEnergy(const double* coord, bool stretched);
// template< typename Impl > void MembraneStretching< Impl, MembraneStretchingAccumulationType::ByTriangle >::computeForces(const double* coord, double* force);
// template< typename Impl > double MembraneStretching< Impl, MembraneStretchingAccumulationType::ByVertex >::computeEnergy(const double* coord, bool stretched);
// template< typename Impl > void MembraneStretching< Impl, MembraneStretchingAccumulationType::ByVertex >::computeForces(const double* coord, double* force);


#endif
