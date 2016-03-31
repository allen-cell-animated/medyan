
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
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

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;


double CylinderExclVolRepulsion::energy(Bead* b1, Bead* b2,
                                        Bead* b3, Bead* b4,
                                        double kRepuls) {
    
    auto c1 = b1->coordinate;
    auto c2 = b2->coordinate;
    auto c3 = b3->coordinate;
    auto c4 = b4->coordinate;
    
    //check if parallel
    if(areParallel(c1, c2, c3, c4)) {
        
        double d = twoPointDistance(c1, c3);
        double invDSquare =  1 / (d * d);
        double energy = kRepuls * invDSquare * invDSquare;
        
        return energy;
    }
    
    //check if in same plane
    if(areInPlane(c1, c2, c3, c4)) {
        
        //slightly move point
        c2 = movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
    }
    
    double a = scalarProduct(c1, c2, c1, c2);
    double b = scalarProduct(c3, c4, c3, c4);
    double c = scalarProduct(c3, c1, c3, c1);
    double d = scalarProduct(c1, c2, c3, c4);
    double e = scalarProduct(c1, c2, c3, c1);
    double f = scalarProduct(c3, c4, c3, c1);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    
    double ATG1 = atan( (a + e)/AA) - atan(e/AA);
    double ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
    double ATG3 = atan((f)/BB) - atan((f - b)/BB);
    double ATG4 = atan((d + f)/FF) - atan((d + f - b)/FF);
    
    double energy = 0.5*kRepuls/JJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
    
    return energy;
}

double CylinderExclVolRepulsion::energy(Bead* b1, Bead* b2,
                                        Bead* b3, Bead* b4,
                                        double kRepuls, double z) {
    
    auto c1 = b1->coordinate;
    auto c2 = b2->coordinate;
    auto c3 = b3->coordinate;
    auto c4 = b4->coordinate;
    
    vector<double> zero (3,0); //Aux zero vector;
    
    vector<double> c1Stretched = {c1[0] + z * b1->force[0],
                                  c1[1] + z * b1->force[1],
                                  c1[2] + z * b1->force[2]};
    
    vector<double> c2Stretched = {c2[0] + z * b2->force[0],
                                  c2[1] + z * b2->force[1],
                                  c2[2] + z * b2->force[2]};
    
    vector<double> c3Stretched = {c3[0] + z * b3->force[0],
                                  c3[1] + z * b3->force[1],
                                  c3[2] + z * b3->force[2]};
    
    vector<double> c4Stretched = {c4[0] + z * b4->force[0],
                                  c4[1] + z * b4->force[1],
                                  c4[2] + z * b4->force[2]};
    
    //check if parallel
    if(areParallel(c1Stretched, c2Stretched,
                   c3Stretched, c4Stretched)) {
        
        double d = twoPointDistance(c1Stretched, c3Stretched);
        double invDSquare =  1 / (d * d);
        double energy =  kRepuls * invDSquare * invDSquare;
        
        return energy;
    }
    
    //check if in same plane
    if(areInPlane(c1Stretched, c2Stretched, c3Stretched, c4Stretched)) {
        //slightly move point
        c2Stretched = movePointOutOfPlane(c1Stretched, c2Stretched,
                                          c3Stretched, c4Stretched, 2, 0.01);
    }
    
    double a = scalarProduct(c1Stretched, c2Stretched, c1Stretched, c2Stretched);
    double b = scalarProduct(c3Stretched, c4Stretched, c3Stretched, c4Stretched);
    double c = scalarProduct(c3Stretched, c1Stretched, c3Stretched, c1Stretched);
    double d = scalarProduct(c1Stretched, c2Stretched, c3Stretched, c4Stretched);
    double e = scalarProduct(c1Stretched, c2Stretched, c3Stretched, c1Stretched);
    double f = scalarProduct(c3Stretched, c4Stretched, c3Stretched, c1Stretched);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    double ATG1 = atan( (a + e)/AA) - atan(e/AA);
    double ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
    double ATG3 = atan((f)/BB) - atan((f - b)/BB);
    double ATG4 = atan((d + f)/FF) - atan((d + f - b)/FF);
    
    double energy = 0.5*kRepuls/JJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
    
    return energy;

}

void CylinderExclVolRepulsion::forces(Bead* b1, Bead* b2,
                                      Bead* b3, Bead* b4,
                                      double kRepuls) {
    
    auto c1 = b1->coordinate;
    auto c2 = b2->coordinate;
    auto c3 = b3->coordinate;
    auto c4 = b4->coordinate;
    
    //check if parallel
    if(areParallel(c1, c2, c3, c4)) {
        
        double d = twoPointDistance(c1, c3);
        double invDSquare =  1/ (d * d);
        double f0 = 4 * kRepuls * invDSquare * invDSquare * invDSquare;
        
        b1->force[0] += - f0 * (c3[0] - c1[0]);
        b1->force[1] += - f0 * (c3[1] - c1[1]);
        b1->force[2] += - f0 * (c3[2] - c1[2]);
        
        b2->force[0] += - f0 * (c4[0] - c2[0]);
        b2->force[1] += - f0 * (c4[1] - c2[1]);
        b2->force[2] += - f0 * (c4[2] - c2[2]);
        
        b3->force[0] += f0 * (c3[0] - c1[0]);
        b3->force[1] += f0 * (c3[1] - c1[1]);
        b3->force[2] += f0 * (c3[2] - c1[2]);
        
        b4->force[0] += f0 * (c4[0] - c2[0]);
        b4->force[1] += f0 * (c4[1] - c2[1]);
        b4->force[2] += f0 * (c4[2] - c2[2]);
        
        return;
    }
    
    //check if in same plane
    if(areInPlane(c1, c2, c3, c4)) {
        
        //slightly move point
        c2 = movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
    }

    double a = scalarProduct(c1, c2, c1, c2);
    double b = scalarProduct(c3, c4, c3, c4);
    double c = scalarProduct(c3, c1, c3, c1);
    double d = scalarProduct(c1, c2, c3, c4);
    double e = scalarProduct(c1, c2, c3, c1);
    double f = scalarProduct(c3, c4, c3, c1);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    double invJJ = 1/JJ;
    
    double ATG1 = atan( (a + e)/AA) - atan(e/AA);
    double ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
    double ATG3 = atan((f)/BB) - atan((f - b)/BB);
    double ATG4 = atan((d + f)/FF) - atan((d + f - b)/FF);
    
    double U = 0.5*kRepuls*invJJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4 );
    
    double A1 = AA*AA/(AA*AA + e*e);
    double A2 = AA*AA/(AA*AA + (a + e)*(a + e));
    
    double E1 = EE*EE/(EE*EE + (a + e - d)*(a + e - d));
    double E2 = EE*EE/(EE*EE + (e - d)*(e - d));
    
    double B1 = BB*BB/(BB*BB + (f - b)*(f - b));;
    double B2 = BB*BB/(BB*BB + f*f);
    
    double F1 = FF*FF/(FF*FF + (d + f - b)*(d + f - b));
    double F2 = FF*FF/(FF*FF + (d + f)*(d + f));
    
    
    double A11 = ATG1/AA;
    double A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) -
                (A2*CC*(a + e))/(AA*AA*AA);
    double A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
    double A14 = (A2*CC)/(AA*AA);
    
    double E11 = ATG2/EE;
    double E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) -
                 (ATG2*GG)/(EE*EE);
    double E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
    double E14 = (E2*GG)/(EE*EE);
    
    double B11 = ATG3/BB;
    double B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*f)/(BB*BB*BB) +
                (B1*DD*(-b + f))/(BB*BB*BB);
    double B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
    double B14 = (B1*DD)/(BB*BB);
    
    double F11 = ATG4/FF;
    double F12 = (F2*(-d - f)*HH)/(FF*FF*FF) + (F1*(-b + d + f)*HH)/(FF*FF*FF) -
                 (ATG4*HH)/(FF*FF);
    double F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
    double F14 = (F1*HH)/(FF*FF);

    b1->force[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (c4[0] - c3[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (c1[0] - c3[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->force[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (c4[1] - c3[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (c1[1] - c3[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->force[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (c4[2] - c3[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (c1[2] - c3[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    b2->force[0] +=  - invJJ*( (c2[0] - c1[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(c4[0] - c3[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(c1[0] - c3[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->force[1] += - invJJ*( (c2[1] - c1[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(c4[1] - c3[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(c1[1] - c3[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->force[2] += - invJJ*( (c2[2] - c1[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(c4[2] - c3[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(c1[2] - c3[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b3->force[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (c4[0] - c3[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[0] - c3[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b3->force[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (c4[1] - c3[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[1] - c3[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
    
    b3->force[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (c4[2] - c3[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[2] - c3[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    b4->force[0] +=  - invJJ*( 0.5*(c2[0] - c1[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (c4[0] - c3[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[0] - c3[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) )  ;
    
    b4->force[1] +=  - invJJ*( 0.5*(c2[1] - c1[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (c4[1] - c3[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[1] - c3[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->force[2] +=  - invJJ*( 0.5*(c2[2] - c1[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (c4[2] - c3[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[2] - c3[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;

}

void CylinderExclVolRepulsion::forcesAux(Bead* b1, Bead* b2,
                                         Bead* b3, Bead* b4,
                                         double kRepuls) {
    
    auto c1 = b1->coordinate;
    auto c2 = b2->coordinate;
    auto c3 = b3->coordinate;
    auto c4 = b4->coordinate;
    
    //check if parallel
    if(areParallel(c1, c2, c3, c4)) {
        
        double d = twoPointDistance(c1, c3);
        double invDSquare =  1/ (d * d);
        double f0 = 4 * kRepuls * invDSquare * invDSquare * invDSquare;
        
        b1->forceAux[0] += - f0 * (c3[0] - c1[0]);
        b1->forceAux[1] += - f0 * (c3[1] - c1[1]);
        b1->forceAux[2] += - f0 * (c3[2] - c1[2]);
        
        b2->forceAux[0] += - f0 * (c4[0] - c2[0]);
        b2->forceAux[1] += - f0 * (c4[1] - c2[1]);
        b2->forceAux[2] += - f0 * (c4[2] - c2[2]);
        
        b3->forceAux[0] += f0 * (c3[0] - c1[0]);
        b3->forceAux[1] += f0 * (c3[1] - c1[1]);
        b3->forceAux[2] += f0 * (c3[2] - c1[2]);
        
        b4->forceAux[0] += f0 * (c4[0] - c2[0]);
        b4->forceAux[1] += f0 * (c4[1] - c2[1]);
        b4->forceAux[2] += f0 * (c4[2] - c2[2]);
        
        return;
    }
    
    //check if in same plane
    if(areInPlane(c1, c2, c3, c4)) {
        
        //slightly move point
        c2 = movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
    }
    
    double a = scalarProduct(c1, c2, c1, c2);
    double b = scalarProduct(c3, c4, c3, c4);
    double c = scalarProduct(c3, c1, c3, c1);
    double d = scalarProduct(c1, c2, c3, c4);
    double e = scalarProduct(c1, c2, c3, c1);
    double f = scalarProduct(c3, c4, c3, c1);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
        
    double invJJ = 1/JJ;
    
    double ATG1 = atan( (a + e)/AA) - atan(e/AA);
    double ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
    double ATG3 = atan((f)/BB) - atan((f - b)/BB);
    double ATG4 = atan((d + f)/FF) - atan((d + f - b)/FF);
    
    double U = 0.5*kRepuls*invJJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4 );
    
    double A1 = AA*AA/(AA*AA + e*e);
    double A2 = AA*AA/(AA*AA + (a + e)*(a + e));
    
    double E1 = EE*EE/(EE*EE + (a + e - d)*(a + e - d));
    double E2 = EE*EE/(EE*EE + (e - d)*(e - d));
    
    double B1 = BB*BB/(BB*BB + (f - b)*(f - b));;
    double B2 = BB*BB/(BB*BB + f*f);
    
    double F1 = FF*FF/(FF*FF + (d + f - b)*(d + f - b));
    double F2 = FF*FF/(FF*FF + (d + f)*(d + f));
    
    
    double A11 = ATG1/AA;
    double A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) -
                 (A2*CC*(a + e))/(AA*AA*AA);
    double A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
    double A14 = (A2*CC)/(AA*AA);
    
    double E11 = ATG2/EE;
    double E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) -
                 (ATG2*GG)/(EE*EE);
    double E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
    double E14 = (E2*GG)/(EE*EE);
    
    double B11 = ATG3/BB;
    double B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*f)/(BB*BB*BB) +
                 (B1*DD*(-b + f))/(BB*BB*BB);
    double B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
    double B14 = (B1*DD)/(BB*BB);
    
    double F11 = ATG4/FF;
    double F12 = (F2*(-d - f)*HH)/(FF*FF*FF) + (F1*(-b + d + f)*HH)/(FF*FF*FF) -
                 (ATG4*HH)/(FF*FF);
    double F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
    double F14 = (F1*HH)/(FF*FF);
    
    b1->forceAux[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (c4[0] - c3[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (c1[0] - c3[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->forceAux[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (c4[1] - c3[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (c1[1] - c3[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->forceAux[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (c4[2] - c3[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (c1[2] - c3[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b2->forceAux[0] +=  - invJJ*( (c2[0] - c1[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(c4[0] - c3[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(c1[0] - c3[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->forceAux[1] += - invJJ*( (c2[1] - c1[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(c4[1] - c3[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(c1[1] - c3[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->forceAux[2] += - invJJ*( (c2[2] - c1[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(c4[2] - c3[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(c1[2] - c3[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b3->forceAux[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (c4[0] - c3[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[0] - c3[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b3->forceAux[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (c4[1] - c3[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[1] - c3[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) / 2;
    
    b3->forceAux[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (c4[2] - c3[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[2] - c3[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b4->forceAux[0] +=  - invJJ*( 0.5*(c2[0] - c1[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (c4[0] - c3[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[0] - c3[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->forceAux[1] +=  - invJJ*( 0.5*(c2[1] - c1[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (c4[1] - c3[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[1] - c3[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->forceAux[2] +=  - invJJ*( 0.5*(c2[2] - c1[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (c4[2] - c3[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[2] - c3[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) );
}

