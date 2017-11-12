
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#include <limits>
typedef std::numeric_limits< double > dbl;

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::vectorize() {

    beadSet = new int[n * Cylinder::getCylinders().size()];
    kstr = new double[Cylinder::getCylinders().size()];
    eql = new double[Cylinder::getCylinders().size()];
    
    int i = 0;
    
    for (auto c: Cylinder::getCylinders()) {
        beadSet[n * i] = c->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = c->getSecondBead()->_dbIndex;

        kstr[i] = c->getMCylinder()->getStretchingConst();
        eql[i] = c->getMCylinder()->getEqLength();
        
        i++;
    }
}

template<class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kstr;
    delete eql;
}


template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    cout.precision(dbl::max_digits10);
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, d);

#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;
    for (auto f: Filament::getFilaments()) {
        
        U_ii = 0;
        
        if (d == 0.0){
            for(auto it : f->getCylinderVector()){
                
                Bead* b1 = it->getFirstBead();
                Bead* b2 = it->getSecondBead();
                double kStretch = it->getMCylinder()->getStretchingConst();
                double eqLength = it->getMCylinder()->getEqLength();
                U_ii += _FFType.energy(b1, b2, kStretch, eqLength);
            }
        }
        else {
            for(auto it : f->getCylinderVector()){
                Bead* b1 = it->getFirstBead();
                Bead* b2 = it->getSecondBead();
                
//                std::cout<<b1->coordinate[0]<<" "<<b1->coordinate[1]<<" "<<b1->coordinate[2]<<" "<<b2->coordinate[0]<<" "<<b2->coordinate[1]<<" "<<b2->coordinate[2]<<" ";
                double kStretch =it->getMCylinder()->getStretchingConst();
                double eqLength = it->getMCylinder()->getEqLength();
                
                U_ii += _FFType.energy(b1, b2, kStretch, eqLength, d);
            }
        }
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_ii != U_ii || U_ii < -1.0) {
            
            //set culprit and return
            _filamentCulprit = f;
            
            U2=-1;
            break;
        }
        else
            U2 += U_ii;
    }
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E S YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }
#endif
    return U_i;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql);
//    for(int i =0;i< Bead::getBeads().size();i++){
//        std::cout<<f[3*i]<<" "<<f[3*i+1]<<" "<<f[3*i+2]<<" ";
//    }
//    std::cout<<endl;

#ifdef CROSSCHECK
    for (auto f: Filament::getFilaments()) {
        
        for(auto it : f->getCylinderVector()){
            
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch =it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();
            
            
            if(cross_checkclass::Aux)
                _FFType.forcesAux(b1, b2, kStretch, eqLength);
            else
                _FFType.forces(b1, b2, kStretch, eqLength);

            
        }
    }
//    for(auto bd:Bead::getBeads())
//        std::cout<<bd->force[0]<<" "<<bd->force[1]<<" "<<bd->force[2]<<" ";
//    std::cout<<endl;

     if(cross_checkclass::Aux)
     {auto state=cross_check::crosscheckAuxforces(f);    std::cout<<"F S YES "<<state<<endl;}
     else
     { auto state=cross_check::crosscheckforces(f);     std::cout<<"F S YES "<<state<<endl;}

#endif

}


///Temlate specializations
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void FilamentStretching<FilamentStretchingHarmonic>::computeForces(double *coord, double *f);
template void FilamentStretching<FilamentStretchingHarmonic>::vectorize();
template void FilamentStretching<FilamentStretchingHarmonic>::deallocate();

