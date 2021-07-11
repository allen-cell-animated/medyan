
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
#include "Mechanics/ForceField/Types.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <unordered_map>

typedef Eigen::Triplet<double> Triplet;

// Forward declarations
class Cylinder;

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
    void vectorizeAllForceFields(const FFCoordinateStartingIndex&);
    /// Deallocation of vectorized memory
    void cleanupAllForceFields();

    /// Compute the energy using all available force fields
    /// @return Returns infinity if there was a problem with a ForceField
    /// energy calculation, such that beads will not be moved to this
    /// problematic configuration.
    /// @param stretched - whether intermediate variables are treated as temporary or not
    template< bool stretched = false >
    
    floatingpoint computeEnergy(floatingpoint *coord, bool verbose = false) const;
    
    
    EnergyReport computeEnergyHRMD(floatingpoint *coord) const;
    
    
    /// Compute the forces of all force fields 
    void computeForces(floatingpoint *coord, std::vector< floatingpoint >& force);
    
    // compute the Hessian matrix if the feature is enabled
    void computeHessian(const std::vector<floatingpoint>& coord, int total_DOF, float delta);
    
    void setCurrBeadMap(const FFCoordinateStartingIndex& si);
    
    // compute the displacement projections along the eigenvectors.
    // Warning: this function only works if all bead coordinates are independent coordinates. Otherwise, out-of-bound access may result.
    void computeProjections(const FFCoordinateStartingIndex&, const std::vector<floatingpoint>& currCoords);
    
    void clearHessian(int a){
        if(a == 0){
            hessianVector.clear();
        }else if(a==1){
            evaluesVector.clear();
            IPRIVector.clear();
            IPRIIVector.clear();
            
        }else{
            projectionsVector.clear();
            tauVector.clear();
        };
    }
    
    vector<floatingpoint> HRMDenergies;
    
    void printculprit();
    
    vector<vector<vector<floatingpoint>>> hessianVector;
    
    vector<Eigen::VectorXcd> evaluesVector;
    vector<Eigen::VectorXcd> IPRIVector;
    vector<Eigen::VectorXcd> IPRIIVector;
    Eigen::VectorXcd evalues;
    Eigen::MatrixXcd evectors;
    vector<floatingpoint> tauVector;
    vector<Eigen::VectorXcd> projectionsVector;
    
    int hessCounter;

    // Map bead pointer to coordinate index in the vectorized data.
    std::unordered_map<Bead*, int> prevBeadMap;
    std::unordered_map<Bead*, int> currBeadMap;
    // Previous coordinates during last eigenvector projection.
    std::vector< floatingpoint > prevCoords;


    vector<string> getinteractionnames(){
        vector<string> temp;
        for (auto &ff : _forceFields)
            for(auto names:ff->getinteractionnames())
            	temp.push_back(names);
        return temp;
    }

#ifdef CUDAACCL
        cudaStream_t  streamF = NULL;
    /// CUDA Copy forces from f to fprev
    void CUDAcopyForces(cudaStream_t  stream, floatingpoint *f, floatingpoint *fprev);
#endif

    /// Compute the load forces on the beads. This does not update the force (xyz) vector
    /// contained by Bead, but updates the loadForce vector which contains precalculated
    /// load values based on the bead's directionality of growth in a filament.
    void computeLoadForces() const;

    // Compute the load forces on the bead for a specific cylinder.
    void computeLoadForce(Cylinder* c, ForceField::LoadForceEnd end) const;
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
