
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_CylinderVolumeFF_h
#define M3SYM_CylinderVolumeFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class CylinderVolumeInteractions;
class Cylinder;

/// An implementation of the ForceField class that calculates
/// Cylinder volume interactions.
class CylinderVolumeFF : public ForceField {
    
private:
    vector <unique_ptr<CylinderVolumeInteractions>>
    _cylinderVolInteractionVector;  ///< Vector of initialized volume interactions
    
    CylinderVolumeInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Initialize the volume forcefields
    CylinderVolumeFF(string& interaction);
    
    virtual string getName() {return "Cylinder Volume";}
    virtual void whoIsCulprit();

    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
