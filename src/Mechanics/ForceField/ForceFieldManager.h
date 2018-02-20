
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
    
    /// Vectorize all interactions involved in calculation
    void vectorizeAllForceFields();
    /// Deallocation of vectorized memory
    void cleanupAllForceFields();
    
    /// Compute the energy using all available force fields
    /// @return Returns infinity if there was a problem with a ForceField
    /// energy calculation, such that beads will not be moved to this
    /// problematic configuration.
    /// @param print - prints detailed info about energies
    double computeEnergy(double *coord, double *f, double d, bool verbose = false);
    
    /// Compute the forces of all force fields 
    void computeForces(double *coord, double *f);
    
    /// Copy forces from f to fprev
    void copyForces(double *f, double *fprev);

#ifdef CUDAACCL
        cudaStream_t  stream;
    /// CUDA Copy forces from f to fprev
    void CUDAcopyForces(cudaStream_t  stream, double *f, double *fprev);
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
#endif
};

#endif
