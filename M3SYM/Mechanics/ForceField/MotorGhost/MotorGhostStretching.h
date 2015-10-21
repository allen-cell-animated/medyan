
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

#ifndef M3SYM_MotorGhostStretching_h
#define M3SYM_MotorGhostStretching_h

#include "common.h"

#include "MotorGhostInteractions.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents a MotorGhost stretching interaction.
template <class MStretchingInteractionType>
class MotorGhostStretching : public MotorGhostInteractions {
    
private:
    MStretchingInteractionType _FFType;

public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Motor Stretching";}
};

#endif
