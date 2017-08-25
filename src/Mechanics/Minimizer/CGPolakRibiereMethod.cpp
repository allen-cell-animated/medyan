

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
    
    startMinimization();
//    long index = 0;
//    long i = 0;
//    for(auto b: Bead::getBeads()) {
//        
//        //flatten indices
//        index = 3 * i;
//        std::cout<<coord[index]<<" "<<coord[index+1]<<" "<<coord[index+2]<<" ";
//        
//        i++;
//    }
//    std::cout<<endl;
    
    FFM.vectorizeAllForceFields();
    
    FFM.computeForces(coord, force);
    FFM.copyForces(forceAux, force);
//    std::cout<<"FORCES"<<endl;
    
//    index=0; i=0;
//    for(auto b: Bead::getBeads()) {
//        
//        //flatten indices
//        index = 3 * i;
//        std::cout<<force[index]<<" "<<force[index+1]<<" "<<force[index+2]<<" ";
//        
//        i++;
//    }
//    std::cout<<endl;
    //compute first gradient
    double curGrad = CGMethod::allFDotF();
  
	int numIter = 0;
//            std::cout<<"CONDITION "<< numIter<<" "<<maxF()<<" "<<GRADTOL<<endl;
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  maxF() > GRADTOL) {
//        std::cout<<"CONDITION "<< numIter<<" "<<maxF()<<" "<<GRADTOL<<endl;
		numIter++;
		double lambda, beta, newGrad, prevGrad;
        
        //find lambda by line search, move beads
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX);

//        std::cout<<endl;
//        std::cout<<"MOVEBEADS"<<endl;
        moveBeads(lambda);
//        index = 0;
//        i = 0;
//        for(auto b: Bead::getBeads()) {
//            
//            //flatten indices
//            index = 3 * i;
//            std::cout<<coord[index]<<" "<<coord[index+1]<<" "<<coord[index+2]<<" ";
//            
//            i++;
//        }
//        std::cout<<endl;
        //compute new forces
        FFM.computeForces(coord, forceAux);
        
        //compute direction
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();
        
        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        
        //shift gradient
        shiftGradient(beta);
        
        //direction reset if not downhill, or no progress made
        if(CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad)) {
            
            shiftGradient(0.0);
            _safeMode = true;
        }
        curGrad = newGrad;
//         index = 0;
//         i = 0;
//        for(auto b: Bead::getBeads()) {
//            
//            //flatten indices
//            index = 3 * i;
//            std::cout<<coord[index]<<" "<<coord[index+1]<<" "<<coord[index+2]<<" ";
//            
//            i++;
//        }

    }
    
    if (numIter >= N) {
        cout << endl;
        
        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;
        
        cout << "Culprit ..." << endl;
        auto b = maxBead();
        if(b != nullptr) b->getParent()->printSelf();
        
        cout << "System energy..." << endl;
        FFM.computeEnergy(coord, force, 0.0, true);
        
        cout << endl;
    }
    
//    cout << "Minimized." << endl;
    
    //final force calculation
    FFM.computeForces(coord, force);
    FFM.copyForces(forceAux, force);
//    index = 0;
//    i = 0;
//    for(auto b: Bead::getBeads()) {
//        
//        //flatten indices
//        index = 3 * i;
//        std::cout<<coord[index]<<" "<<coord[index+1]<<" "<<coord[index+2]<<" ";
//        
//        i++;
//    }
//    std::cout<<endl;
    endMinimization();
    FFM.computeLoadForces();

    
    FFM.cleanupAllForceFields();
}
