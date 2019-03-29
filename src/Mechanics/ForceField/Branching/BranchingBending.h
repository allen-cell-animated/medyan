
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BranchingBending_h
#define MEDYAN_BranchingBending_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction maintaining a BranchingPoint angle (~70 for Arp2/3)
template <class BBendingInteractionType>
class BranchingBending : public BranchingInteractions {
    
private:
    BBendingInteractionType _FFType;
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    double *kbend;
    double *eqt;
#ifdef CUDAACCL
    int * gpu_beadSet;
    double * gpu_kbend;
    double *gpu_eqt;
    int * gpu_params;
    CUDAvars cvars;
    double *F_i;
#endif
public:
    
    ///Array describing indexed set of interactions
    ///this is a 4-bead potential
    const static int n = 4;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual const string getName() {return "Branching Bending";}
};

#endif
