
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

#include "CGSteepestDescent.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Mechanics/Minimizer/CGDataCopy.hpp"
#include "Output.h"
#include "Structure/Bead.h"

MinimizationResult SteepestDescent::minimize(
    ForceFieldManager &FFM, floatingpoint GRADTOL,
    floatingpoint MAXDIST, floatingpoint LAMBDAMAX,
    floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
    string _LINESEARCHALGORITHM,
    bool steplimit
) {
    
    MinimizationResult result;

        //number of steps
        int N;
        if (steplimit) {
            int beadMaxStep = 5 * Bead::numBeads();
            N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
        } else {
            N = numeric_limits<int>::max();
        }

    startMinimization(GRADTOL);
    FFM.vectorizeAllForceFields(initCGMethodData(*this, GRADTOL));

    FFM.computeForces(coord.data(), force);
    searchDir = force;
    auto maxForce = maxF();
    bool isForceBelowTol = forceBelowTolerance();

    result.energiesBefore = FFM.computeEnergyHRMD(coord.data());

        int numIter = 0;
        while (/* Iteration criterion */  numIter < N &&
               /* Gradient tolerance  */  ! isForceBelowTol) {

            numIter++;
            floatingpoint lambda;

            //find lambda by line search, move beads
            bool *dummy = nullptr;
            lambda = backtrackingLineSearch(FFM, MAXDIST, maxForce, LAMBDAMAX,
                    LAMBDARUNNINGAVERAGEPROBABILITY, dummy, dummy);
            moveAlongSearchDir(lambda);

        //compute new forces
        FFM.computeForces(coord.data(), force);
        maxForce = maxF();
        isForceBelowTol = forceBelowTolerance();

        //shift gradient
        shiftSearchDir(0.0);
    }

        if (numIter >= N) {
            cout << endl;

            cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
            cout << "Maximum force in system = " << maxF() << endl;

            cout << "Culprit ..." << endl;
            auto b = maxBead();
            if (b != nullptr) b->getParent()->printSelf();

            cout << "System energy..." << endl;
            FFM.computeEnergy(coord.data(), true);

            cout << endl;
        }

    result.energiesAfter = FFM.computeEnergyHRMD(coord.data());

    //final force calculation
    FFM.computeForces(coord.data(), force);
    searchDir = force;
    FFM.computeLoadForces();

    // Copy the coordinate and force data back to the system
    copyFromCGMethodData(*this);
    endMinimization();

    FFM.cleanupAllForceFields();

    return result;
}

