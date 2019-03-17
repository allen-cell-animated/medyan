
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

#include "CGSteepestDescent.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"

    void SteepestDescent::minimize(ForceFieldManager &FFM, double GRADTOL,
                                   double MAXDIST, double LAMBDAMAX, bool steplimit) {
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

        FFM.computeForces(coord, force);
        Bead::getDbData().forcesAux = Bead::getDbData().forces;

        int numIter = 0;
        while (/* Iteration criterion */  numIter < N &&
                                          /* Gradient tolerance  */  maxF() > GRADTOL) {

            numIter++;
            double lambda;

            //find lambda by line search, move beads
            bool *dummy = nullptr;
            lambda = backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
            moveBeads(lambda);

            //compute new forces
            FFM.computeForces(coord, forceAux);

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
            FFM.computeEnergy(coord, force, 0.0, true);

            cout << endl;
        }

        //final force calculation
        FFM.computeForces(coord, force);
        Bead::getDbData().forcesAux = Bead::getDbData().forces;
        FFM.computeLoadForces();
        endMinimization();

        FFM.cleanupAllForceFields();
    }

