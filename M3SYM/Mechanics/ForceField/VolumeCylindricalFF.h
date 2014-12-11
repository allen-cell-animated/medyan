
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_VolumeCylindricalFF_h
#define M3SYM_VolumeCylindricalFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class CylinderVolumeInteractions;

/// An implementation of the ForceField class that calculates Cylinder
/// volume interactions.
class VolumeCylindricalFF : public ForceField {
    
private:
    vector <unique_ptr<CylinderVolumeInteractions>> _cylinderVolInteractionVector;  ///< Vector of initialized volume interactions
    
public:
    /// Initialize the volume forcefields
    VolumeCylindricalFF(string& interaction);
    
    virtual string getName() {return "Volume cylindrical";}
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif
