
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

#ifndef MEDYAN_MotorGhostStretching_h
#define MEDYAN_MotorGhostStretching_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

#include "MotorGhostInteractions.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents a MotorGhost stretching interaction.
template <class MStretchingInteractionType>
class MotorGhostStretching : public MotorGhostInteractions {
    
private:
    MStretchingInteractionType _FFType;

    int *beadSet;
    
    ///Array describing the constants in calculation
    double *kstr;
    double *eql;
    double *pos1;
    double *pos2;
    double *stretchforce;

#ifdef CUDAACCL
    int * gpu_beadSet;
    double * gpu_kstr;
    double *gpu_eql;
    int * gpu_params;
    double *gpu_pos1;
    double *gpu_pos2;
    double *F_i;
    double *gpu_Mstretchforce;
    cudaStream_t stream = NULL;
#endif
    
public:
    ///Array describing indexed set of interactions
    ///For linkers, this is a 4-bead potential
    const static int n = 4;
    
    ///< Constructor
    MotorGhostStretching () {}
    ~MotorGhostStretching () {}
    
    virtual void vectorize();
    virtual void deallocate();
    
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);


    virtual const string getName() {return "MotorGhost Stretching";}

    virtual void assignforcemags();

};

#endif
