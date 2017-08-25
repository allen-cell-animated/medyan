
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

#include "MotorGhostStretchingHarmonic.h"
#include "MotorGhostStretching.h"
#include "MotorGhost.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MotorGhostStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2) {
    
    int n = MotorGhostStretching<MotorGhostStretchingHarmonic>::n;
    int nint = MotorGhost::getMotorGhosts().size();
    
    double *coord1, *coord2, *coord3, *coord4, dist, U_i;
    double *v1 = new double[3];
    double *v2 = new double[3];
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];
        
        midPointCoordinate(v1, coord1, coord2, pos1[i]);
        midPointCoordinate(v2, coord3, coord4, pos2[i]);
        
        dist = twoPointDistance(v1, v2) - eql[i];
        U_i = 0.5 * kstr[i] * dist * dist;
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    delete v1;
    delete v2;
    
    return U;
}

double MotorGhostStretchingHarmonic::energy(double *coord, double * f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2, double d){
    
    int n = MotorGhostStretching<MotorGhostStretchingHarmonic>::n;
    int nint = MotorGhost::getMotorGhosts().size();
    
    double *coord1, *coord2, *coord3, *coord4, *f1, *f2, *f3, *f4, dist, U_i;
    double *v1 = new double[3];
    double *v2 = new double[3];
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];
        
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        f4 = &f[3 * beadSet[n * i + 3]];
        
        midPointCoordinateStretched(v1, coord1, f1, coord2, f2, pos1[i], d);
        midPointCoordinateStretched(v2, coord3, f3, coord4, f4, pos2[i], d);
        
        dist = twoPointDistance(v1,  v2) - eql[i];
        U_i = 0.5 * kstr[i] * dist * dist;
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    delete v1;
    delete v2;
    
    return U;
    
}
void MotorGhostStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql, double *pos1, double *pos2){
    
    
    int n = MotorGhostStretching<MotorGhostStretchingHarmonic>::n;
    int nint = MotorGhost::getMotorGhosts().size();
    
    double *coord1, *coord2, *coord3, *coord4, dist, invL;
    double *v1 = new double[3];
    double *v2 = new double[3];
    
    double f0, *f1, *f2, *f3, *f4;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];
        
        midPointCoordinate(v1, coord1, coord2, pos1[i]);
        midPointCoordinate(v2, coord3, coord4, pos2[i]);
        
        dist = twoPointDistance(v1, v2) ;
        invL = 1 / dist;
        
        f0 = kstr[i] * ( dist - eql[i] ) * invL;
        
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        f4 = &f[3 * beadSet[n * i + 3]];
        
        //force on i
        f1[0] +=   -f0 * ( v1[0] - v2[0] ) * (1 - pos1[i]);
        f1[1] +=   -f0 * ( v1[1] - v2[1] ) * (1 - pos1[i]);
        f1[2] +=   -f0 * ( v1[2] - v2[2] ) * (1 - pos1[i]);
        
        // force i+1
        f2[0] +=   -f0 * ( v1[0] - v2[0] ) * (pos1[i]);
        f2[1] +=   -f0 * ( v1[1] - v2[1] ) * (pos1[i]);
        f2[2] +=   -f0 * ( v1[2] - v2[2] ) * (pos1[i]);
        
        //force on j
        f3[0] +=   f0 * ( v1[0] - v2[0] ) * (1 - pos2[i]);
        f3[1] +=   f0 * ( v1[1] - v2[1] ) * (1 - pos2[i]);
        f3[2] +=   f0 * ( v1[2] - v2[2] ) * (1 - pos2[i]);
        
        // force j+1
        f4[0] +=   f0 * ( v1[0] - v2[0] ) * (pos2[i]);
        f4[1] +=   f0 * ( v1[1] - v2[1] ) * (pos2[i]);
        f4[2] +=   f0 * ( v1[2] - v2[2] ) * (pos2[i]);
//             std::cout<<"MOTOR "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
    }
    delete v1;
    delete v2;
}
