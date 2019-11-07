#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolume_hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolume_hpp

#include <memory> // unique_ptr
#include <string>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/NeighborListImpl.h"
#include "Structure/Bead.h"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

//FORWARD DECLARATIONS
class Cylinder;

/// Represents an excuded volume interaction between a triangle and a cylinder (bead).
template < typename InteractionType >
class TriangleBeadExclVolume : public ForceField {
    
private:
    InteractionType _FFType;
    std::unique_ptr<TriangleFilBeadNL> _neighborList;  ///< Neighbor list of triangle-bead

    // Culprit elements
    const Triangle* triangleCulprit_ = nullptr;
    const Bead*     beadCulprit_     = nullptr;

public:
    ///Constructor
    TriangleBeadExclVolume() :
        _neighborList(std::make_unique< TriangleFilBeadNL >(
            SysParams::Mechanics().MemBeadVolumeCutoff,
            SysParams::Mechanics().MemBeadVolumeCutoffMech
        ))
    { }

    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override;
    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override;

    virtual void computeLoadForces() override;
    virtual void computeLoadForce(Cylinder* c, LoadForceEnd end) const override;

    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {
        return { _neighborList.get() };
    }

    virtual std::string getName() override {return "Triangle Bead Excluded Volume";}

    virtual void whoIsCulprit() override {
        LOG(INFO) << "Printing the culprit triangle and bead...";
        triangleCulprit_->printSelf();
        beadCulprit_    ->printSelf();
    }

    // Useless overrides
    virtual void vectorize() override {}
    virtual void cleanup() override {}
};

#endif
