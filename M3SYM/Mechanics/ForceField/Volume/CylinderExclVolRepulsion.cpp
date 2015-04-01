
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
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

using namespace mathfunc;


double CylinderExclVolRepulsion::energy(Bead* b1, Bead* b2,
                                        Bead* b3, Bead* b4,
                                        double kRepuls) {
    
    if(ifParallel(b1->coordinate, b2->coordinate,
                  b3->coordinate, b4->coordinate)) {
        
        double d = twoPointDistance(b1->coordinate, b3->coordinate);
        double invDSquare =  1 / (d * d);
        double energy = kRepuls * invDSquare * invDSquare;
        
        return energy;
    }
    
    double a = scalarProduct(b1->coordinate, b2->coordinate,
                             b1->coordinate, b2->coordinate);
    double b = scalarProduct(b3->coordinate, b4->coordinate,
                             b3->coordinate, b4->coordinate);
    double c = scalarProduct(b3->coordinate, b1->coordinate,
                             b3->coordinate, b1->coordinate);
    double d = scalarProduct(b1->coordinate, b2->coordinate,
                             b3->coordinate, b4->coordinate);
    double e = scalarProduct(b1->coordinate, b2->coordinate,
                             b3->coordinate, b1->coordinate);
    double f = scalarProduct(b3->coordinate, b4->coordinate,
                             b3->coordinate, b1->coordinate);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    if (abs(JJ) < 1e-10 || JJ != JJ) {

        auto v = movePointOutOfPlane(b1->coordinate, b2->coordinate,
                                     b3->coordinate, b4->coordinate, 4, 1.0);
        
        a = scalarProduct(v, b2->coordinate, v, b2->coordinate);
        b = scalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate);
        c = scalarProduct(b3->coordinate, v, b3->coordinate, v);
        d = scalarProduct(v, b2->coordinate, b3->coordinate, b4->coordinate);
        e = scalarProduct(v, b2->coordinate, b3->coordinate, v);
        f = scalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, v);
        
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - f*f);
        
        CC = d*e - a*f;
        DD = b*e - d*f;
        
        EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
        
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - f*CC;
        
    }
    
    double ATG1 = atan( (a + e)/AA) - atan(e/AA);
    double ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
    double ATG3 = atan((f)/BB) - atan((f - b)/BB);
    double ATG4 = atan((d + f)/FF) - atan((d + f - b)/FF);
    
    double energy = 0.5*kRepuls/JJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
    return energy;
}

double CylinderExclVolRepulsion::energy(Bead* b1, Bead* b2,
                                        Bead* b3, Bead* b4,
                                        double kRepuls, double lambda) {
    
    if(ifParallel(b1->coordinate, b2->coordinate,
                  b3->coordinate, b4->coordinate)) {
        
        double d = twoPointDistanceStretched(b1->coordinate, b1->force,
                                             b3->coordinate, b3->force, lambda);
        double invDSquare =  1/ (d * d);
        double energy =  kRepuls * invDSquare * invDSquare;
        
        return energy;
    }
    
    double a = scalarProductStretched(b1->coordinate, b1->force,
                                      b2->coordinate, b2->force,
                                      b1->coordinate, b1->force,
                                      b2->coordinate, b2->force, lambda);
    double b = scalarProductStretched(b3->coordinate, b3->force,
                                      b4->coordinate, b4->force,
                                      b3->coordinate, b3->force,
                                      b4->coordinate, b4->force, lambda);
    double c = scalarProductStretched(b3->coordinate, b3->force,
                                      b1->coordinate, b1->force,
                                      b3->coordinate, b3->force,
                                      b1->coordinate, b1->force, lambda);
    double d = scalarProductStretched(b1->coordinate, b1->force,
                                      b2->coordinate, b2->force,
                                      b3->coordinate, b3->force,
                                      b4->coordinate, b4->force, lambda);
    double e = scalarProductStretched(b1->coordinate, b1->force,
                                      b2->coordinate, b2->force,
                                      b3->coordinate, b3->force,
                                      b1->coordinate, b1->force, lambda);
    double f = scalarProductStretched(b3->coordinate, b3->force,
                                      b4->coordinate, b4->force,
                                      b3->coordinate, b3->force,
                                      b1->coordinate, b1->force, lambda);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    if (abs(JJ) < 1e-10 || JJ != JJ){
        
        auto v = movePointOutOfPlane(b1->coordinate, b2->coordinate,
                                     b3->coordinate, b4->coordinate, 4, 1.0);
        
        a = scalarProductStretched(v, b1->force, b2->coordinate, b2->force,
                                   v, b1->force, b2->coordinate, b2->force, lambda);
        b = scalarProductStretched(b3->coordinate, b3->force, b4->coordinate, b4->force,
                                   b3->coordinate, b3->force, b4->coordinate, b4->force, lambda);
        c = scalarProductStretched(b3->coordinate, b3->force, v, b1->force,
                                   b3->coordinate, b3->force, v, b1->force, lambda);
        d = scalarProductStretched(v, b1->force, b2->coordinate, b2->force,
                                   b3->coordinate, b3->force, b4->coordinate, b4->force, lambda);
        e = scalarProductStretched(v, b1->force, b2->coordinate, b2->force,
                                   b3->coordinate, b3->force, v, b1->force, lambda);
        f = scalarProductStretched(b3->coordinate, b3->force, b4->coordinate, b4->force,
                                   b3->coordinate, b3->force, v, b1->force, lambda);
        
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - f*f);
        
        CC = d*e - a*f;
        DD = b*e - d*f;
        
        EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
        
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - f*CC;
    }
    
    
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
    
    if(ifParallel(b1->coordinate, b2->coordinate, b3->coordinate, b4->coordinate)) {
        
        double d = twoPointDistance(b1->coordinate, b3->coordinate);
        double invDSquare =  1/ (d * d);
        double f0 = 4 * kRepuls * invDSquare * invDSquare * invDSquare;
        
        b1->force[0] += - f0 * (b3->coordinate[0] - b1->coordinate[0]);
        b1->force[1] += - f0 * (b3->coordinate[1] - b1->coordinate[1]);
        b1->force[2] += - f0 * (b3->coordinate[2] - b1->coordinate[2]);
        
        b2->force[0] += - f0 * (b4->coordinate[0] - b2->coordinate[0]);
        b2->force[1] += - f0 * (b4->coordinate[1] - b2->coordinate[1]);
        b2->force[2] += - f0 * (b4->coordinate[2] - b2->coordinate[2]);
        
        b3->force[0] += f0 * (b3->coordinate[0] - b1->coordinate[0]);
        b3->force[1] += f0 * (b3->coordinate[1] - b1->coordinate[1]);
        b3->force[2] += f0 * (b3->coordinate[2] - b1->coordinate[2]);
        
        b4->force[0] += f0 * (b4->coordinate[0] - b2->coordinate[0]);
        b4->force[1] += f0 * (b4->coordinate[1] - b2->coordinate[1]);
        b4->force[2] += f0 * (b4->coordinate[2] - b2->coordinate[2]);
        
        return;
    }

    double a = scalarProduct(b1->coordinate, b2->coordinate,
                             b1->coordinate, b2->coordinate);
    double b = scalarProduct(b3->coordinate, b4->coordinate,
                             b3->coordinate, b4->coordinate);
    double c = scalarProduct(b3->coordinate, b1->coordinate,
                             b3->coordinate, b1->coordinate);
    double d = scalarProduct(b1->coordinate, b2->coordinate,
                             b3->coordinate, b4->coordinate);
    double e = scalarProduct(b1->coordinate, b2->coordinate,
                             b3->coordinate, b1->coordinate);
    double f = scalarProduct(b3->coordinate, b4->coordinate,
                             b3->coordinate, b1->coordinate);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    if (abs(JJ) < 1e-10 || JJ != JJ){
        
        auto v = movePointOutOfPlane(b1->coordinate, b2->coordinate,
                                     b3->coordinate, b4->coordinate, 4, 1.0);
        
        a = scalarProduct(v, b2->coordinate, v, b2->coordinate);
        b = scalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate);
        c = scalarProduct(b3->coordinate, v, b3->coordinate, v);
        d = scalarProduct(v, b2->coordinate, b3->coordinate, b4->coordinate);
        e = scalarProduct(v, b2->coordinate, b3->coordinate, v);
        f = scalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, v);
        
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - f*f);
        
        CC = d*e - a*f;
        DD = b*e - d*f;
        
        EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
        
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - f*CC;
    }
    
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

    b1->force[0] +=  - 0.5*invJJ*( (b2->coordinate[0] - b1->coordinate[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinate[0] - b3->coordinate[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinate[0] - b3->coordinate[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->force[1] +=  - 0.5*invJJ*( (b2->coordinate[1] - b1->coordinate[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinate[1] - b3->coordinate[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinate[1] - b3->coordinate[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->force[2] +=  - 0.5*invJJ*( (b2->coordinate[2] - b1->coordinate[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinate[2] - b3->coordinate[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinate[2] - b3->coordinate[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    b2->force[0] +=  - invJJ*( (b2->coordinate[0] - b1->coordinate[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinate[0] - b3->coordinate[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinate[0] - b3->coordinate[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->force[1] += - invJJ*( (b2->coordinate[1] - b1->coordinate[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinate[1] - b3->coordinate[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinate[1] - b3->coordinate[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->force[2] += - invJJ*( (b2->coordinate[2] - b1->coordinate[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinate[2] - b3->coordinate[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinate[2] - b3->coordinate[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b3->force[0] +=  - 0.5*invJJ*( (b2->coordinate[0] - b1->coordinate[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinate[0] - b3->coordinate[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinate[0] - b3->coordinate[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b3->force[1] +=  - 0.5*invJJ*( (b2->coordinate[1] - b1->coordinate[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinate[1] - b3->coordinate[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinate[1] - b3->coordinate[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
    
    b3->force[2] +=  - 0.5*invJJ*( (b2->coordinate[2] - b1->coordinate[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinate[2] - b3->coordinate[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinate[2] - b3->coordinate[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    b4->force[0] +=  - invJJ*( 0.5*(b2->coordinate[0] - b1->coordinate[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinate[0] - b3->coordinate[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinate[0] - b3->coordinate[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) )  ;
    
    b4->force[1] +=  - invJJ*( 0.5*(b2->coordinate[1] - b1->coordinate[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinate[1] - b3->coordinate[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinate[1] - b3->coordinate[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->force[2] +=  - invJJ*( 0.5*(b2->coordinate[2] - b1->coordinate[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinate[2] - b3->coordinate[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinate[2] - b3->coordinate[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;

}

void CylinderExclVolRepulsion::forcesAux(Bead* b1, Bead* b2,
                                         Bead* b3, Bead* b4,
                                         double kRepuls) {
    
    if(ifParallel(b1->coordinateAux, b2->coordinateAux,
                  b3->coordinateAux, b4->coordinateAux)) {
        
        double d = twoPointDistance(b1->coordinateAux, b3->coordinateAux);
        double invDSquare =  1/ (d * d);
        double f0 = 4 * kRepuls * invDSquare * invDSquare * invDSquare;
        
        b1->forceAux[0] += - f0 * (b3->coordinateAux[0] - b1->coordinateAux[0]);
        b1->forceAux[1] += - f0 * (b3->coordinateAux[1] - b1->coordinateAux[1]);
        b1->forceAux[2] += - f0 * (b3->coordinateAux[2] - b1->coordinateAux[2]);
        
        b2->forceAux[0] += - f0 * (b4->coordinateAux[0] - b2->coordinateAux[0]);
        b2->forceAux[1] += - f0 * (b4->coordinateAux[1] - b2->coordinateAux[1]);
        b2->forceAux[2] += - f0 * (b4->coordinateAux[2] - b2->coordinateAux[2]);
        
        b3->forceAux[0] += f0 * (b3->coordinateAux[0] - b1->coordinateAux[0]);
        b3->forceAux[1] += f0 * (b3->coordinateAux[1] - b1->coordinateAux[1]);
        b3->forceAux[2] += f0 * (b3->coordinateAux[2] - b1->coordinateAux[2]);
        
        b4->forceAux[0] += f0 * (b4->coordinateAux[0] - b2->coordinateAux[0]);
        b4->forceAux[1] += f0 * (b4->coordinateAux[1] - b2->coordinateAux[1]);
        b4->forceAux[2] += f0 * (b4->coordinateAux[2] - b2->coordinateAux[2]);
        
        return;
    }
    
    double a = scalarProduct(b1->coordinateAux, b2->coordinateAux,
                             b1->coordinateAux, b2->coordinateAux);
    double b = scalarProduct(b3->coordinateAux, b4->coordinateAux,
                             b3->coordinateAux, b4->coordinateAux);
    double c = scalarProduct(b3->coordinateAux, b1->coordinateAux,
                             b3->coordinateAux, b1->coordinateAux);
    double d = scalarProduct(b1->coordinateAux, b2->coordinateAux,
                             b3->coordinateAux, b4->coordinateAux);
    double e = scalarProduct(b1->coordinateAux, b2->coordinateAux,
                             b3->coordinateAux, b1->coordinateAux);
    double f = scalarProduct(b3->coordinateAux, b4->coordinateAux,
                             b3->coordinateAux, b1->coordinateAux);
    
    double AA = sqrt(a*c - e*e);
    double BB = sqrt(b*c - f*f);
    
    double CC = d*e - a*f;
    double DD = b*e - d*f;
    
    double EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
    double FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
    
    
    double GG = d*d - a*b - CC;
    double HH = CC + GG - DD;
    double JJ = c*(GG + CC) + e*DD - f*CC;
    
    if (abs(JJ) < 1e-10 || JJ != JJ){
        
        auto v = movePointOutOfPlane(b1->coordinateAux, b2->coordinateAux,
                                     b3->coordinateAux, b4->coordinateAux, 4, 1.0);
        
        a = scalarProduct(v, b2->coordinateAux, v, b2->coordinateAux);
        b = scalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, b4->coordinateAux);
        c = scalarProduct(b3->coordinateAux, v, b3->coordinateAux, v);
        d = scalarProduct(v, b2->coordinateAux, b3->coordinateAux, b4->coordinateAux);
        e = scalarProduct(v, b2->coordinateAux, b3->coordinateAux, v);
        f = scalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, v);
        
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - f*f);
        
        CC = d*e - a*f;
        DD = b*e - d*f;
        
        EE = sqrt( a*(b + c - 2*f) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + f)*(d + f) );
        
        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - f*CC;
    }
    
    
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
    
    b1->forceAux[0] +=  - 0.5*invJJ*( (b2->coordinateAux[0] - b1->coordinateAux[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinateAux[0] - b3->coordinateAux[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinateAux[0] - b3->coordinateAux[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->forceAux[1] +=  - 0.5*invJJ*( (b2->coordinateAux[1] - b1->coordinateAux[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinateAux[1] - b3->coordinateAux[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinateAux[1] - b3->coordinateAux[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b1->forceAux[2] +=  - 0.5*invJJ*( (b2->coordinateAux[2] - b1->coordinateAux[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinateAux[2] - b3->coordinateAux[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinateAux[2] - b3->coordinateAux[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b2->forceAux[0] +=  - invJJ*( (b2->coordinateAux[0] - b1->coordinateAux[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinateAux[0] - b3->coordinateAux[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinateAux[0] - b3->coordinateAux[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->forceAux[1] += - invJJ*( (b2->coordinateAux[1] - b1->coordinateAux[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinateAux[1] - b3->coordinateAux[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinateAux[1] - b3->coordinateAux[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->forceAux[2] += - invJJ*( (b2->coordinateAux[2] - b1->coordinateAux[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinateAux[2] - b3->coordinateAux[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinateAux[2] - b3->coordinateAux[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b3->forceAux[0] +=  - 0.5*invJJ*( (b2->coordinateAux[0] - b1->coordinateAux[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinateAux[0] - b3->coordinateAux[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinateAux[0] - b3->coordinateAux[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b3->forceAux[1] +=  - 0.5*invJJ*( (b2->coordinateAux[1] - b1->coordinateAux[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinateAux[1] - b3->coordinateAux[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinateAux[1] - b3->coordinateAux[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) / 2;
    
    b3->forceAux[2] +=  - 0.5*invJJ*( (b2->coordinateAux[2] - b1->coordinateAux[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinateAux[2] - b3->coordinateAux[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinateAux[2] - b3->coordinateAux[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b4->forceAux[0] +=  - invJJ*( 0.5*(b2->coordinateAux[0] - b1->coordinateAux[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinateAux[0] - b3->coordinateAux[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinateAux[0] - b3->coordinateAux[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->forceAux[1] +=  - invJJ*( 0.5*(b2->coordinateAux[1] - b1->coordinateAux[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinateAux[1] - b3->coordinateAux[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinateAux[1] - b3->coordinateAux[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->forceAux[2] +=  - invJJ*( 0.5*(b2->coordinateAux[2] - b1->coordinateAux[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinateAux[2] - b3->coordinateAux[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinateAux[2] - b3->coordinateAux[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) );
}

