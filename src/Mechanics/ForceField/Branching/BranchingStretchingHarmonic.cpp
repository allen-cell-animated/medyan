
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

#include "BranchingStretchingHarmonic.h"
#include "BranchingStretching.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BranchingStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                           double *kstr, double *eql, double *pos){
    
    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, dist, U_i;
    double *v1 = new double[3];
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        
        midPointCoordinate(v1, coord1, coord2, pos[i]);
        dist = twoPointDistance(v1, coord3) - eql[i];
        
        U_i = 0.5 * kstr[i] * dist * dist;
        
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    delete v1;
    return U;
}

double BranchingStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                           double *kstr, double *eql, double *pos, double d){
    
    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *f1, *f2, *f3, dist, U_i;
    double *v1 = new double[3];
    double *vzero = new double[3]; vzero[0] = 0; vzero[1] = 0; vzero[2] = 0;

    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        
        
        midPointCoordinateStretched(v1, coord1, f1, coord2, f2, pos[i], d);
        dist = twoPointDistanceStretched(v1, vzero, coord3, f3, d) - eql[i];
        
        U_i = 0.5 * kstr[i] * dist * dist;
        
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    delete v1; vzero;
    return U;
}

void BranchingStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                         double *kstr, double *eql, double *pos){
    
    
    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *f1, *f2, *f3, dist, invL, f0;
    double *v1 = new double[3];

    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        
        
        midPointCoordinate(v1, coord1, coord2, pos[i]);
        dist = twoPointDistance(v1, coord3) - eql[i];
        
        invL = 1 / dist;
        f0 = kstr[i] * ( dist - eql[i]) * invL;
    
        f1[0] +=  -f0 * ( coord3[0] - v1[0] ) * (pos[i] - 1);
        f1[1] +=  -f0 * ( coord3[1] - v1[1] ) * (pos[i] - 1);
        f1[2] +=  -f0 * ( coord3[2] - v1[2] ) * (pos[i] - 1);
        
        // force i+1
        f2[0] +=  f0 * ( coord3[0] - v1[0] ) * pos[i];
        f2[1] +=  f0 * ( coord3[1] - v1[1] ) * pos[i];
        f2[2] +=  f0 * ( coord3[2] - v1[2] ) * pos[i];
        
        //force on j
        f3[0] +=  -f0 * ( coord3[0] - v1[0] );
        f3[1] +=  -f0 * ( coord3[1] - v1[1] );
        f3[2] +=  -f0 * ( coord3[2] - v1[2] );
        
    }
    delete v1;
}
