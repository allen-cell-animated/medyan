
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

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"

void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL,
                            double MAXDIST, double LAMBDAMAX, bool steplimit){
    
    //number of steps
    int N;
    if(steplimit) {
        int beadMaxStep = 5 * Bead::numBeads();
        N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
    }
    else {
        N = numeric_limits<int>::max();
    }
    
	FFM.computeForces();
    startMinimization();

    //compute first gradient
    double curGrad = CGMethod::allFDotF();

    std::cout<<"printing beads & forces"<<endl;
    for(auto b:Bead::getBeads()){
        auto x = b->coordinate;
        std::cout<<b->getID()<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<b->force[0]<<" "
                ""<<b->force[1]<<" "<<b->force[2]<<endl;
    }
    std::cout<<"printed beads & forces"<<endl;

	int numIter = 0;
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  maxF() > GRADTOL) {

		numIter++;
		double lambda, beta, newGrad, prevGrad;
        std::cout<<"SL maxF "<<maxF()<<endl;
        
        //find lambda by line search, move beads
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX);
        
        moveBeads(lambda); setBeads();
        
        //compute new forces
        FFM.computeForcesAux();
        std::cout<<"MB printing beads & forces L "<<lambda<<endl;
        for(auto b:Bead::getBeads()){
            auto x = b->coordinate;
            std::cout<<b->getID()<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<b->forceAux[0]<<" "
                    ""<<b->forceAux[1]<<" "<<b->forceAux[2]<<endl;
        }
        std::cout<<"MB printed beads & forces"<<endl;
        //compute direction
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();
        
        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        
        //update prev forces
        FFM.computeForcesAuxP();
        
        //shift gradient
        shiftGradient(beta);
        std::cout<<"Shift Gradient "<<beta<<endl;

        //direction reset if not downhill, or no progress made
        if(CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad)) {
            
            shiftGradient(0.0);
            _safeMode = true;
            std::cout<<"Shift Gradient 0.0"<<endl;
        }
        curGrad = newGrad;
        std::cout<<"Maximum Force"<<maxF()<<endl;
    }
    std::cout<<"Total iterations "<<numIter<<endl;

    if (numIter >= N) {
        cout << endl;
        
        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;
        
        cout << "Culprit ..." << endl;
        auto b = maxBead();
        if(b != nullptr) b->getParent()->printSelf();
        
        cout << "System energy..." << endl;
        FFM.computeEnergy(0.0, true);
        
        cout << endl;
    }
    
    //final force calculation
    FFM.computeForces();
    std::cout<<"End Minimization......"<<endl;
    FFM.computeLoadForces();
    std::cout<<"printing beads & forces"<<endl;
    for(auto b:Bead::getBeads()){
        auto x = b->coordinate;
        std::cout<<b->getID()<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<b->force[0]<<" "
                ""<<b->force[1]<<" "<<b->force[2]<<endl;
    }
    std::cout<<"printed beads & forces"<<endl;
}
