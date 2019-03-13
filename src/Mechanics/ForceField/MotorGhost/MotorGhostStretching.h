
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

#ifndef MEDYAN_MotorGhostStretching_h
#define MEDYAN_MotorGhostStretching_h

#include "common.h"
#include "CUDAcommon.h"

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
    floatingpoint *kstr;
    floatingpoint *eql;
    floatingpoint *pos1;
    floatingpoint *pos2;
    floatingpoint *stretchforce;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kstr;
    floatingpoint *gpu_eql;
    int * gpu_params;
    floatingpoint *gpu_pos1;
    floatingpoint *gpu_pos2;
    floatingpoint *F_i;
    floatingpoint *gpu_Mstretchforce;
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
    
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);


    virtual const string getName() {return "MotorGhost Stretching";}

    virtual void assignforcemags();

};

#endif
