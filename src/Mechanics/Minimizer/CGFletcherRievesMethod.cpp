
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

#include "CGFletcherRievesMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"
#include "Structure/Bead.h"
    void FletcherRieves::minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
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

        //compute first gradient
        floatingpoint curGrad = CGMethod::allFDotF();

        int numIter = 0;
        while (/* Iteration criterion */  numIter < N &&
                                          /* Gradient tolerance  */  maxF() > GRADTOL) {
            numIter++;
            floatingpoint lambda, beta, newGrad;

            //temporary
            bool *dummy = nullptr;
            //find lambda by line search, move beads
            lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
                               : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
            moveBeads(lambda);

            //compute new forces
            FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data());

            //compute direction
            newGrad = CGMethod::allFADotFA();

            //Fletcher-Rieves update
            beta = newGrad / curGrad;

            //shift gradient
            shiftGradient(beta);

            //direction reset if not downhill or no progress made
            if (CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad)) {
                shiftGradient(0.0);
                _safeMode = true;
            }

            curGrad = newGrad;
        }

        if (numIter >= N) {
            cout << endl;

            cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
            cout << "Maximum force in system = " << maxF() << endl;

            cout << "Culprit ..." << endl;
            auto b = maxBead();
            if (b != nullptr) b->getParent()->printSelf();

            cout << "System energy..." << endl;
            FFM.computeEnergy(Bead::getDbData().coords.data(), Bead::getDbData().forces.data(), 0.0, true);

            cout << endl;
        }

        //final force calculation
        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
        Bead::getDbData().forcesAux = Bead::getDbData().forces;
        FFM.computeLoadForces();
        endMinimization();

        FFM.cleanupAllForceFields();
    }

