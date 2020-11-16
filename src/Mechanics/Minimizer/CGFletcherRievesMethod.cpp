
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

#include "CGFletcherRievesMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"
#include "Structure/Bead.h"
MinimizationResult FletcherRieves::minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                                  floatingpoint MAXDIST, floatingpoint LAMBDAMAX,
                                  floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
                                  string _LINESEARCHALGORITHM,
                                  bool steplimit) {

    MinimizationResult result;

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
        result.energiesBefore = FFM.computeEnergyHRMD(Bead::getDbData().coords.data());

        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
        Bead::getDbData().forcesAux = Bead::getDbData().forces;
        auto maxForce = maxF();

        //compute first gradient
        floatingpoint curGrad = CGMethod::allFDotF();

        int numIter = 0;
        while (/* Iteration criterion */  numIter < N &&
               /* Gradient tolerance  */  maxForce > GRADTOL) {
            numIter++;
            floatingpoint lambda, beta, newGrad;

            //temporary
            bool *dummy = nullptr;
            //find lambda by line search, move beads
            lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, maxForce,
            		LAMBDAMAX, dummy, dummy)
                               : backtrackingLineSearch(FFM, MAXDIST, maxForce, LAMBDAMAX,
                                       LAMBDARUNNINGAVERAGEPROBABILITY, dummy, dummy);
            moveBeads(lambda);

            //compute new forces
            FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data());
            maxForce = maxF();

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
            FFM.computeEnergy(Bead::getDbData().coords.data(), true);

            cout << endl;
        }

        result.energiesAfter = FFM.computeEnergyHRMD(Bead::getDbData().coords.data());

        //final force calculation
        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
        Bead::getDbData().forcesAux = Bead::getDbData().forces;
        FFM.computeLoadForces();
        endMinimization();

        FFM.cleanupAllForceFields();

    return result;
}

