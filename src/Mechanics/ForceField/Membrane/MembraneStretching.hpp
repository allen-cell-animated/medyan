#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "utility.h"

enum class MembraneStretchingAccumulationType {
    ByTriangle, ByVertex
};

template< typename Impl, MembraneStretchingAccumulationType >
class MembraneStretching: public ForceField {
private:
    Impl _impl;

    // Culprit membrane
    const Membrane* membraneCulprit_ = nullptr;

    // Force buffer
    std::vector< floatingpoint > forceBuffer_;

public:
    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override;
    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override;

    virtual string getName() override { return "Membrane Stretching"; }

    virtual void whoIsCulprit() override {
        if(membraneCulprit_) {
            LOG(INFO) << "Printing culprit membrane...";
            membraneCulprit_->printSelf();
        } else {
            LOG(ERROR) << "Membrane culprit is not set.";
        }
    }

    const auto& getForceBuffer() const { return forceBuffer_; }

    // Useless overrides
    virtual void vectorize() override {}
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return std::vector<NeighborList*>(); }
};

#endif
