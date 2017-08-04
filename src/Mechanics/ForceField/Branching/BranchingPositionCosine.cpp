
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

#include <cmath>
#include <math.h>

#include "BranchingPositionCosine.h"
#include "BranchingPosition.h"

#include "BranchingPoint.h"

#include "Bead.h"

#include "MathFunctions.h"


using namespace mathfunc;

double BranchingPositionCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kpos, double *pos){
    
    
    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = n * BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, X, D, XD, xd, theta, eqTheta, dTheta, U_i;
    double *mp = new double[3];
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        
        midPointCoordinate(mp, coord1, coord2, pos[i]);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));
        
        XD = X * D;
        
        xd = sqrt(scalarProduct(mp, coord2, mp, coord3));
        
        theta = safeacos(xd / XD);
        eqTheta = 0.5*M_PI;
        dTheta = theta-eqTheta;
        
        U_i = kpos[i] * ( 1 - cos(dTheta) );
        
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    delete mp;
    return U;
}

double BranchingPositionCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kpos, double *pos, double d){
    
    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = n * BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *f1, *f2, *f3, X, D, XD, xd, theta, eqTheta, dTheta, U_i;
    double *mp = new double[3];
    double *vzero = new double[3]; vzero[0] = 0; vzero[1] = 0; vzero[2] = 0;
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
    

        midPointCoordinateStretched(mp, coord1, f1, coord2, f2, pos[i], d);
        X = sqrt(scalarProductStretched(mp, vzero, coord2, f2, mp, vzero, coord2, f2, d));
        D = sqrt(scalarProductStretched(mp, vzero, coord3, f3, mp, vzero, coord3, f3, d));
        
        XD = X * D;
        
        xd = sqrt(scalarProductStretched(mp, vzero, coord2, f2, mp, vzero, coord3, f3, d));
        
        theta = safeacos(xd / XD);
        eqTheta = 0.5*M_PI;
        dTheta = theta-eqTheta;
        
        U_i = kpos[i] * ( 1 - cos(dTheta) );
        
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    delete mp;
    delete vzero;
    return U;
}

void BranchingPositionCosine::forces(double *coord, double *f, int *beadSet,
                                     double *kpos, double *pos){
    
    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = n * BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *f1, *f2, *f3, X, D, XD, xd, invX, invD, position, A, B, C, k, theta, eqTheta, dTheta, U_i;
    double *mp = new double[3];
    
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        
        midPointCoordinate(mp, coord1, coord2, pos[i]);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));
        
        XD = X * D;
        
        xd = sqrt(scalarProduct(mp, coord2, mp, coord3));
        
        invX = 1/X;
        invD = 1/D;
        A = invX*invD;
        B = invX*invX;
        C = invD*invD;
        
        theta = safeacos(xd / XD);
        eqTheta = 0.5*M_PI;
        dTheta = theta-eqTheta;
        
        position = pos[i];
        
        k =  kpos[i] * A * sin(dTheta)/sin(theta);
        
        //bead 1
        f1[0] +=  k * (1-position)* (- (1-position)*(coord2[0] - coord1[0]) - (coord3[0] - (1-position)*coord1[0] - position*coord2[0])
                            + xd *(B*(1-position)*(coord2[0] - coord1[0]) + C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0])) );
        
        f1[1] +=  k * (1-position)* (- (1-position)*(coord2[1] - coord1[1]) - (coord3[1] - (1-position)*coord1[1] - position*coord2[1])
                            + xd *(B*(1-position)*(coord2[1] - coord1[1]) + C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1])) );
        
        f1[2] +=  k * (1-position)* (- (1-position)*(coord2[2] - coord1[2]) - (coord3[2] - (1-position)*coord1[2] - position*coord2[2])
                            + xd *(B*(1-position)*(coord2[2] - coord1[2]) + C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2])) );
        
        //bead 2
        
        f2[0] +=  k * (- position*(1-position)*(coord2[0] - coord1[0]) + (1-position)*(coord3[0]- (1-position)*coord1[0] - position*coord2[0])
            + xd *( (1-position)*B*(1-position)*(coord2[0] - coord1[0]) - position*C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0])) );
        
        f2[1] +=  k * (- position*(1-position)*(coord2[1] - coord1[1]) + (1-position)*(coord3[1]- (1-position)*coord1[1] - position*coord2[1])
            + xd *( (1-position)*B*(1-position)*(coord2[1] - coord1[1]) - position*C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1])) );
        
        f2[2] +=  k * (- position*(1-position)*(coord2[2] - coord1[2]) + (1-position)*(coord3[2]- (1-position)*coord1[2] - position*coord2[2])
            + xd *( (1-position)*B*(1-position)*(coord2[2] - coord1[2]) - position*C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2])) );
        
        //bead3
        
        f3[0] +=  k * ( (1-position)*(coord2[0] - coord1[0]) - xd * C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) );
        f3[1] +=  k * ( (1-position)*(coord2[1] - coord1[1]) - xd * C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) );
        f3[2] +=  k * ( (1-position)*(coord2[2] - coord1[2]) - xd * C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) );
        
    }
    delete mp;
}
