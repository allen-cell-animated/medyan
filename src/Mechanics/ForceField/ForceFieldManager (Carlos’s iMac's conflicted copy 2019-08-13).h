
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

#ifndef MEDYAN_ForceFieldManager_h
#define MEDYAN_ForceFieldManager_h

#include <vector>

#include "common.h"

#include "ForceField.h"
#include "Bead.h"

/// A class to store and iterate over all [ForceFields](@ref ForceField).
/*!
 *  The ForceFieldManager is used to store all [ForceFields](@ref ForceField)
 *  initialized by the system, as well as iterate over these potentials and calculate
 *  total forces and energies. This class contains functions for the said calculations.
 */
class ForceFieldManager {

friend class CGMethod;

public:
     vector<ForceField*> _forceFields; ///< All forcefields in the system

     static ForceField* _culpritForceField;

    /// Vectorize all interactions involved in calculation
    void vectorizeAllForceFields();
    /// Deallocation of vectorized memory
    void cleanupAllForceFields();

    /// Compute the energy using all available force fields
    /// @return Returns infinity if there was a problem with a ForceField
    /// energy calculation, such that beads will not be moved to this
    /// problematic configuration.
    /// @param print - prints detailed info about energies
    floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d, bool verbose = false, bool HRMDbool = false);

    /// Compute the forces of all force fields
    void computeForces(floatingpoint *coord, floatingpoint *f);
    
    // compute the Hessian matrix if the feature is enabled
    void computeHessian(floatingpoint *coord, floatingpoint *f, int total_DOF, float delta);
    
    void clearHessian(){
        hessianVector.clear();
        tauVector.clear();
    }
    
    vector<floatingpoint> HRMDenergies;
    
    void printculprit(floatingpoint* force);
    
    vector<vector<vector<floatingpoint>>> hessianVector;
    
    vector<floatingpoint> tauVector;

#ifdef CUDAACCL
        cudaStream_t  streamF = NULL;
    /// CUDA Copy forces from f to fprev
    void CUDAcopyForces(cudaStream_t  stream, floatingpoint *f, floatingpoint *fprev);
#endif

    /// Compute the load forces on the beads. This does not update the force (xyz) vector
    /// contained by Bead, but updates the loadForce vector which contains precalculated
    /// load values based on the bead's directionality of growth in a filament.
    void computeLoadForces();
#ifdef CROSSCHECK
    /// Reset the forces of all objects
    void resetForces();
#endif
#ifdef CUDAACCL
    vector<int> blocksnthreads;
    int *gpu_nint;
    //@{
    vector<int> bntaddvec2;
    int *gpu_params;
    vector<int> params;
    //@}

#endif
    void assignallforcemags();

private:
    chrono::high_resolution_clock::time_point tbegin, tend;
};

#endif
