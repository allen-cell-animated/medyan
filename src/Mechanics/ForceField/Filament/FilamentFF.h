
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

#ifndef MEDYAN_FilamentFF_h
#define MEDYAN_FilamentFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class FilamentInteractions;
class Filament;

/// An implementation of the ForceField class that calculates Filament interactions.
class FilamentFF : public ForceField {
 
private:
    vector<unique_ptr<FilamentInteractions>>
    _filamentInteractionVector; ///< Vector of initialized filament interactions
    
protected:
    FilamentInteractions* _culpritInteraction; ///< Culprit in case of error
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    FilamentFF(string& stretching, string& bending, string& twisting);
    
    virtual void vectorize();
    virtual void cleanup();
    
    virtual string getName() {return "Filament";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
};

#endif
