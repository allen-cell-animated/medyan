
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CGMethod.h"

#include "ForceFieldManager.h"
#include "Bead.h"

#include "BoundaryElement.h"

#include "BoundaryElementImpl.h" //added by jl135 to accomodate CG for spherical



double CGMethod::allFDotF()
{
    double g = 0;
    for(auto b: Bead::getBeads())
        g += b->FDotF();

    //added by jl135, no need to dot product because there is no normal vector involved
    for(auto be: BoundaryElement::getBoundaryElements()){
    	g+= be->boundary_force * be->boundary_force;
    }

    return g;
}

double CGMethod::allFADotFA()
{
    double g = 0;
	for(auto b: Bead::getBeads())
        g += b->FADotFA();
    
	//added by jl135, no need to dot product because there is no normal vector involved
	for(auto be: BoundaryElement::getBoundaryElements()){
		g+= be->boundary_forceAux * be->boundary_forceAux;
	}

    return g;
}

double CGMethod::allFADotFAP()
{
    double g = 0;
    for(auto b: Bead::getBeads())
        g += b->FADotFAP();
    
    //added by jl135, no need to dot product because there is no normal vector involved

    for(auto be: BoundaryElement::getBoundaryElements()){
    	g+= be->boundary_forceAux * be->boundary_forceAuxP;

    }


    return g;
}

double CGMethod::allFDotFA()
{
    double g = 0;
	for(auto b: Bead::getBeads())
        g += b->FDotFA();
    
	//added by jl135, no need to dot product because there is no normal vector involved
	for(auto be: BoundaryElement::getBoundaryElements()){
		g+= be->boundary_force * be->boundary_forceAux;
	}

    return g;
}

double CGMethod::maxF() {
    
    double maxF = 0;
    

    //calc max force
    for(auto b: Bead::getBeads())
        maxF = max(maxF, sqrt(b->FADotFA()));

    return maxF;
}

//added by jl135, maxFb is to calculate the maximum boundary tension
double CGMethod::maxFb() {

    double maxFb = 0;

	//calc max boundary tension
	for(auto be: BoundaryElement::getBoundaryElements())

		maxFb = max(maxFb, abs(be->boundary_force));
		//abs is needed since boundary_tension will be negative after passing down the equilibrium radius (shrinking)


	    return maxFb;
}

Bead* CGMethod::maxBead() {
    
    double maxF = 0;
    double force;
    Bead* retbead = nullptr;
    
    //calc max force
    for(auto b: Bead::getBeads()) {
        
        force = b->FADotFA();
        
        if(force >= maxF) {
            retbead = b;
            maxF = force;
        }
    }
    return retbead;
}



void CGMethod::startMinimization() {

    for(auto b: Bead::getBeads())
        b->coordinateB = b->coordinateP = b->coordinate;

    //added by jl135, to minimize the total force on boundary
    for(auto be: BoundaryElement::getBoundaryElements()) {
    	SphereBoundaryElement* bs = (SphereBoundaryElement*)be;
    	bs->_radiusB = bs->_radiusP = bs->_radius;
    }
}



void CGMethod::moveBeads(double d)
{
	for(auto b: Bead::getBeads()) {
        
        b->coordinate[0] = b->coordinate[0] + d * b->force[0];
        b->coordinate[1] = b->coordinate[1] + d * b->force[1];
        b->coordinate[2] = b->coordinate[2] + d * b->force[2];
	}

	//added by jl135, to minimize the total force on boundary
    for(auto be: BoundaryElement::getBoundaryElements()){

    	SphereBoundaryElement* bs = (SphereBoundaryElement*)be;
    	bs->_radius = bs->_radius + d* be->boundary_force; ///
    }
}

void CGMethod::resetBeads() {
    
    for(auto b: Bead::getBeads())
        b->coordinate = b->coordinateP;

    //added by jl135, to minimize the total force on boundary
    for(auto be: BoundaryElement::getBoundaryElements()){

    	SphereBoundaryElement* bs = (SphereBoundaryElement*)be;
    	bs->_radius = bs->_radiusP;
    }
}
void CGMethod::setBeads() {
    
    for(auto b: Bead::getBeads())
        b->coordinateP = b->coordinate;

    //added by jl135, to minimize the total force on boundary
    for(auto be: BoundaryElement::getBoundaryElements()){

        SphereBoundaryElement* bs = (SphereBoundaryElement*)be;
        bs->_radiusP = bs->_radius;

    }
}

void CGMethod::shiftGradient(double d)
{
	for(auto b: Bead::getBeads()) {
        
        b->force[0] = b->forceAux[0] + d * b->force[0];
        b->force[1] = b->forceAux[1] + d * b->force[1];
        b->force[2] = b->forceAux[2] + d * b->force[2];
	}

    //added by jl135, to minimize the total force on boundary
    for(auto be: BoundaryElement::getBoundaryElements())
    	be->boundary_force = be->boundary_forceAux + d * be->boundary_force;

}

void CGMethod::printForces()
{
	cout << "Print Forces" << endl;
    for(auto b: Bead::getBeads()) {
        
		for (int i = 0; i<3; i++)
            cout << b->coordinate[i] << "  "<<
                    b->force[i] <<"  "<<b->forceAux[i]<<endl;
	}
    cout << "End of Print Forces" << endl;

    cout << "Print Boundary Forces" << endl;
    for(auto be: BoundaryElement::getBoundaryElements())
    	cout << ((SphereBoundaryElement*) be)->_radius << "   "<< be->boundary_force << "   "<<be->boundary_forceAux <<endl;
    cout << "End of Print Boundary Forces" <<endl;
}

double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                                double LAMBDAMAX) {

    double f = maxF();
    
    //return zero if no forces
    if(f == 0.0) return 0.0;
    
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / f); //CG method at least the line search, depends on bead not boundary tension,
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFA();
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || lambda <= LAMBDATOL)
            return 0.0;
    }
}

double CGMethod::safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                                    double LAMBDAMAX) {
    
    //reset safe mode
    _safeMode = false;
    
    //calculate first lambda
    double lambda = LAMBDAMAX;
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= 0.0) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        //just shake if we cant find an energy min,
        //so we dont get stuck
        if(lambda <= 0.0 || lambda <= LAMBDATOL)
            return MAXDIST / maxF();
    }


}
