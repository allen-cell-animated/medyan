
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
#include "Output.h"
#include "Structure/Bead.h"

    void SteepestDescent::minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                                   floatingpoint MAXDIST, floatingpoint LAMBDAMAX, bool steplimit) {
        //number of steps
        int N;
        if (steplimit) {
            int beadMaxStep = 5 * Bead::numBeads();
            N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
        } else {
            N = numeric_limits<int>::max();
        }

        startMinimization();
        FFM.vectorizeAllForceFields();

        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
        Bead::getDbData().forcesAux = Bead::getDbData().forces;

        int numIter = 0;
        while (/* Iteration criterion */  numIter < N &&
                                          /* Gradient tolerance  */  maxF() > GRADTOL) {

            numIter++;
            floatingpoint lambda;

            //find lambda by line search, move beads
            bool *dummy = nullptr;
            lambda = backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
            moveBeads(lambda);

            //compute new forces
            FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data());

            //shift gradient
            shiftGradient(0.0);
        }

        if (numIter >= N) {
            cout << endl;

            cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
            cout << "Maximum force in system = " << maxF() << endl;

            cout << "Culprit ..." << endl;
            auto b = maxBead();
            if (b != nullptr) b->getParent()->printSelf();

            cout << "System energy..." << endl;
            FFM.computeEnergy(Bead::getDbData().coords.data(), true);

            cout << endl;
        }

        //final force calculation
        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
        Bead::getDbData().forcesAux = Bead::getDbData().forces;
        FFM.computeLoadForces();
        endMinimization();

        FFM.cleanupAllForceFields();
    }

