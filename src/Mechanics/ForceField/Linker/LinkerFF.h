
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_LinkerFF_h
#define MEDYAN_LinkerFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class LinkerInteractions;
class Linker;

/// An implementation of the ForceField class that calculates Linker
/// stretching, bending, and twisting.
class LinkerFF : public ForceField {
    
private:
    vector<unique_ptr<LinkerInteractions>>
    _linkerInteractionVector; ///< Vector of initialized linker interactions
    
protected:
    /// The culprit in the case of an error
    static LinkerInteractions* _culpritInteraction;
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    LinkerFF(string& stretching, string& bending, string& twisting );
    
    virtual string getName() {return "Linker";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
    
};

#endif
