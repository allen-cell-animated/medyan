
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

/*
 Cylinder 1: x1------x2
 Cylinder 2: y1------y2
 
 
 x = x2-x1
 y = y2-y1
 z = x1-y1
 
 a = x.x;
 b = y.y;
 c = z.z;
 d = x.y;
 e = x.z;
 f = y.z;
 
 */

#include "CylinderExclVolRepulsion.h"
#include "CylinderExclVolume.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;


double CylinderExclVolRepulsion::energy(double *coord, double *force, int *beadSet, double *krep) {
    
    
    double *c1, *c2, *c3, *c4, d, invDSquare, U, U_i;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
    double ATG1, ATG2, ATG3, ATG4;
    
    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;
    
    U_i = 0;
    
    for (int i = 0; i < nint; nint++) {
        
        c1 = &coord[3 * beadSet[n * i]];
        c2 = &coord[3 * beadSet[n * i + 1]];
        c3 = &coord[3 * beadSet[n * i + 2]];
        c4 = &coord[3 * beadSet[n * i + 3]];
        
        
        //check if parallel
        if(areParallel(c1, c2, c3, c4)) {
    
            d = twoPointDistance(c1, c3);
            invDSquare =  1 / (d * d);
            U_i = krep[i] * invDSquare * invDSquare;
    
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return TODO
                return -1;
            }
            continue;
        }
    
        //check if in same plane
        if(areInPlane(c1, c2, c3, c4)) {
    
            //slightly move point
            movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
        }
    
        a = scalarProduct(c1, c2, c1, c2);
        b = scalarProduct(c3, c4, c3, c4);
        c = scalarProduct(c3, c1, c3, c1);
        d = scalarProduct(c1, c2, c3, c4);
        e = scalarProduct(c1, c2, c3, c1);
        F = scalarProduct(c3, c4, c3, c1);
    
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - F*F);
    
        CC = d*e - a*F;
        DD = b*e - d*F;
    
        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );
    
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - F*CC;
    
    
        ATG1 = atan( (a + e)/AA) - atan(e/AA);
        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
        ATG3 = atan((F)/BB) - atan((F - b)/BB);
        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);
        
        U_i = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return TODO
            
            return -1;
        }
        U += U_i;
    }
    return U;
}

double CylinderExclVolRepulsion::energy(double *coord, double *f, int *beadSet, double *krep, double z) {

    
    double d, invDSquare, U, U_i, *f1, *f2, *f3, *f4;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
    double ATG1, ATG2, ATG3, ATG4;
    
    double *c1 = new double[3];
    double *c2 = new double[3];
    double *c3 = new double[3];
    double *c4 = new double[3];
    
    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;
    
    U_i = 0;
    
    for (int i = 0; i < nint; nint++) {
        
        memcpy(c1, &coord[3 * beadSet[n * i]], 3 * sizeof(double));
        memcpy(c2, &coord[3 * beadSet[n * i + 1]], 3 * sizeof(double));
        memcpy(c3, &coord[3 * beadSet[n * i + 2]], 3 * sizeof(double));
        memcpy(c4, &coord[3 * beadSet[n * i + 3]], 3 * sizeof(double));
        
        //stretch coords
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        f4 = &f[3 * beadSet[n * i + 3]];
        
        c1[0] = c1[0] + z * f1[0];
        c1[1] = c1[1] + z * f1[1];
        c1[2] = c1[2] + z * f1[2];
        c2[0] = c2[0] + z * f2[0];
        c2[1] = c2[1] + z * f2[1];
        c2[2] = c2[2] + z * f2[2];
        c3[0] = c3[0] + z * f3[0];
        c3[1] = c3[1] + z * f3[1];
        c3[2] = c3[2] + z * f3[2];
        c4[0] = c4[0] + z * f4[0];
        c4[1] = c4[1] + z * f4[1];
        c4[2] = c4[2] + z * f4[2];
        
        //check if parallel
        if(areParallel(c1, c2, c3, c4)) {
            
            d = twoPointDistance(c1, c3);
            invDSquare =  1 / (d * d);
            U_i = krep[i] * invDSquare * invDSquare;
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return TODO
                return -1;
            }
            continue;
        }
        
        //check if in same plane
        if(areInPlane(c1, c2, c3, c4)) {
            
            //slightly move point
            movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
        }
        
        a = scalarProduct(c1, c2, c1, c2);
        b = scalarProduct(c3, c4, c3, c4);
        c = scalarProduct(c3, c1, c3, c1);
        d = scalarProduct(c1, c2, c3, c4);
        e = scalarProduct(c1, c2, c3, c1);
        F = scalarProduct(c3, c4, c3, c1);
        
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - F*F);
        
        CC = d*e - a*F;
        DD = b*e - d*F;
        
        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );
        
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - F*CC;
        
        
        ATG1 = atan( (a + e)/AA) - atan(e/AA);
        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
        ATG3 = atan((F)/BB) - atan((F - b)/BB);
        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);
        
        U_i = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return TODO
            
            return -1;
        }
        U += U_i;
    }
    return U;
}

void CylinderExclVolRepulsion::forces(double *coord, double *f, int *beadSet, double *krep) {
   
    
    double *c1, *c2, *c3, *c4, d, invDSquare, U, *f1, *f2, *f3, *f4;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ, invJJ;
    double ATG1, ATG2, ATG3, ATG4;
    double A1, A2, E1, E2, B1, B2, F1, F2, A11, A12, A13, A14;
    double E11, E12, E13, E14, B11, B12, B13, B14, F11, F12, F13, F14;
    
    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;
    
    for (int i = 0; i < nint; nint++) {
        
        c1 = &coord[3 * beadSet[n * i]];
        c2 = &coord[3 * beadSet[n * i + 1]];
        c3 = &coord[3 * beadSet[n * i + 2]];
        c4 = &coord[3 * beadSet[n * i + 3]];
        
        //stretch coords
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        f4 = &f[3 * beadSet[n * i + 3]];
        
        //check if parallel
        if(areParallel(c1, c2, c3, c4)) {
            
            d = twoPointDistance(c1, c3);
            invDSquare =  1 / (d * d);
            U = krep[i] * invDSquare * invDSquare;
            
            continue;
        }
        
        //check if in same plane
        if(areInPlane(c1, c2, c3, c4)) {
            
            //slightly move point
            movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
        }
        
        a = scalarProduct(c1, c2, c1, c2);
        b = scalarProduct(c3, c4, c3, c4);
        c = scalarProduct(c3, c1, c3, c1);
        d = scalarProduct(c1, c2, c3, c4);
        e = scalarProduct(c1, c2, c3, c1);
        F = scalarProduct(c3, c4, c3, c1);
        
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - F*F);
        
        CC = d*e - a*F;
        DD = b*e - d*F;
        
        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );
        
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - F*CC;
        
        invJJ = 1/JJ;
        
        ATG1 = atan( (a + e)/AA) - atan(e/AA);
        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
        ATG3 = atan((F)/BB) - atan((F - b)/BB);
        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);
        
        U = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
        
        A1 = AA*AA/(AA*AA + e*e);
        A2 = AA*AA/(AA*AA + (a + e)*(a + e));
        
        E1 = EE*EE/(EE*EE + (a + e - d)*(a + e - d));
        E2 = EE*EE/(EE*EE + (e - d)*(e - d));
        
        B1 = BB*BB/(BB*BB + (F - b)*(F - b));;
        B2 = BB*BB/(BB*BB + F*F);
        
        F1 = FF*FF/(FF*FF + (d + F - b)*(d + F - b));
        F2 = FF*FF/(FF*FF + (d + F)*(d + F));
        
        A11 = ATG1/AA;
        A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) -
                (A2*CC*(a + e))/(AA*AA*AA);
        A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
        A14 = (A2*CC)/(AA*AA);
        
        E11 = ATG2/EE;
        E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) -
              (ATG2*GG)/(EE*EE);
        E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
        E14 = (E2*GG)/(EE*EE);
        
        B11 = ATG3/BB;
        B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*F)/(BB*BB*BB) +
                (B1*DD*(-b + F))/(BB*BB*BB);
        B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
        B14 = (B1*DD)/(BB*BB);
        
        F11 = ATG4/FF;
        F12 = (F2*(-d - F)*HH)/(FF*FF*FF) + (F1*(-b + d + F)*HH)/(FF*FF*FF) -
               (ATG4*HH)/(FF*FF);
        F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
        F14 = (F1*HH)/(FF*FF);
    
        f1[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[0] - c3[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[0] - c3[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
        f1[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[1] - c3[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[1] - c3[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
        f1[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[2] - c3[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[2] - c3[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
        f2[0] +=  - invJJ*( (c2[0] - c1[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[0] - c3[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[0] - c3[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );
    
        f2[1] += - invJJ*( (c2[1] - c1[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[1] - c3[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[1] - c3[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );
    
        f2[2] += - invJJ*( (c2[2] - c1[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[2] - c3[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[2] - c3[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );
    
        f3[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[0] - c3[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[0] - c3[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
        f3[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[1] - c3[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[1] - c3[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
    
        f3[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[2] - c3[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[2] - c3[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
        f4[0] +=  - invJJ*( 0.5*(c2[0] - c1[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[0] - c3[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[0] - c3[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) )  ;
    
        f4[1] +=  - invJJ*( 0.5*(c2[1] - c1[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[1] - c3[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[1] - c3[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;
    
        f4[2] +=  - invJJ*( 0.5*(c2[2] - c1[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[2] - c3[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[2] - c3[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;
    }
}
