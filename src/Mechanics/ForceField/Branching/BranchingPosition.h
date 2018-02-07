
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

#ifndef MEDYAN_BranchingPosition_h
#define MEDYAN_BranchingPosition_h

#include "common.h"
#include "CUDAcommon.h"
#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction fixing a Cylinder anchored by a BranchingPoint on the parent.
template <class BStretchingInteractionType>
class BranchingPosition : public BranchingInteractions {
    
private:
    BStretchingInteractionType _FFType;
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    double *kpos;
    double *pos;
#ifdef CUDAACCL
    int * gpu_beadSet;
    double * gpu_kpos;
    double *gpu_pos;
    int * gpu_params;
    CUDAvars cvars;
    double *F_i;
#endif
    
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 3-bead potential
    const static int n = 3;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual const string getName() {return "Branching Position";}
};

#endif

