
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
    LinkerInteractions* _culpritInteraction;
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    LinkerFF(string& stretching, string& bending, string& twisting );
    
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    virtual void cleanup();

    virtual string getName() {return "Linker";}
    virtual void whoIsCulprit();
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, bool stretched = false) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
    //Assigns stretchforces for ratechangeimpl
    virtual void assignforcemags();
    
};

#endif
