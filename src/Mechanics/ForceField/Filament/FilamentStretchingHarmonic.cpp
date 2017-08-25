
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

#include "FilamentStretchingHarmonic.h"
#include "FilamentStretching.h"

#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double FilamentStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql){

    
    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();
    
    double *coord1, *coord2, *coord3, *coord4, dist, U_i;
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        
        dist = twoPointDistance(coord1, coord2) - eql[i];
        U_i = 0.5 * kstr[i] * dist * dist;
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            FilamentInteractions::_filamentCulprit = (Filament*)(Cylinder::getCylinders()[i]->getParent());
            
            return -1;
        }
        
        U += U_i;
    }
    
    return U;
}

double FilamentStretchingHarmonic::energy(double *coord, double * f, int *beadSet,
                                          double *kstr, double *eql, double d){
    
    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();
    
    double *coord1, *coord2, *coord3, *coord4, *f1, *f2, *f3, *f4, dist;
    double *v1 = new double[3];
    double *v2 = new double[3];
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        
        dist = twoPointDistanceStretched(coord1, f1,  coord2, f2, d) - eql[i];
        U += 0.5 * kstr[i] * dist * dist;
    }
    delete v1;
    delete v2;
    
    return U;
}

void FilamentStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                        double *kstr, double *eql){
    
    
    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();
    
    double *coord1, *coord2, *coord3, *coord4, dist, invL;
    double f0, *f1, *f2, *f3, *f4;
//    auto xxx=Cylinder::getCylinders();
    for(int i = 0; i < nint; i += 1) {
//        std::cout<<xxx.size()<<" "<<beadSet[n*i]<<" "<<beadSet[n*i+1]<<endl;
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        
        dist = twoPointDistance(coord1, coord2);
        invL = 1 / dist;
        
        f0 = kstr[i] * ( dist - eql[i] ) * invL;
        
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        
        f1[0] +=  f0 * ( coord1[0] - coord2[0] );
        f1[1] +=  f0 * ( coord1[1] - coord2[1] );
        f1[2] +=  f0 * ( coord1[2] - coord2[2] );
        
        // force i-1
        f2[0] +=  f0 * ( coord2[0] - coord1[0] );
        f2[1] +=  f0 * ( coord2[1] - coord1[1] );
        f2[2] +=  f0 * ( coord2[2] - coord1[2] );
//        std::cout<<"STRETCHING "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
    }
}


