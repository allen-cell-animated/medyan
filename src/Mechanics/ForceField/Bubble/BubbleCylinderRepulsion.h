
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

#ifndef MEDYAN_BubbleCylinderRepulsion_h
#define MEDYAN_BubbleCylinderRepulsion_h

#include <vector>

#include "common.h"

#include "BubbleInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive interaction between a Bubble and Cylinder.
template <class BRepulsionInteractionType>
class BubbleCylinderRepulsion : public BubbleInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    BubbleCylinderNL* _neighborList; ///<Neighbor list of Bubble-Cylinder
    
    int *beadSet;
    int *nneighbors;
    int *bubbleSet;
//    int *nintvec;
    
    ///Array describing the constants in calculation
    double *krep;
    double *slen;
    double *radius;
    int nint = 0;
public:
    
    ///Array describing indexed set of interactions
    ///For bubble, this is a 1-bead potential + 1 fixed bubble bead
    const static int n = 1;
    
    /// Constructor
    BubbleCylinderRepulsion() {
        _neighborList = new BubbleCylinderNL(SysParams::Mechanics().BubbleCutoff);
    }
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    //virtual void computeForcesAux(double *coord, double *f);
    
    virtual void computeLoadForces();
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Bubble-Cylinder Repulsion";}
};

#endif

