//
//  CylindricalVolRepulsion.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
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


double CylinderExclVolRepulsion::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double kRepuls)
{
    double a = ScalarProduct(b1->coordinate, b2->coordinate, b1->coordinate, b2->coordinate);
    double b = ScalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate);
    double c = ScalarProduct(b3->coordinate, b1->coordinate, b3->coordinate, b1->coordinate);
    double d = ScalarProduct(b1->coordinate, b2->coordinate, b3->coordinate, b4->coordinate);
    double e = ScalarProduct(b1->coordinate, b2->coordinate, b3->coordinate, b1->coordinate);
    double f = ScalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b1->coordinate);
    
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
    
    return 0.5*kRepuls/JJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4 );
    
    
    
}

double CylinderExclVolRepulsion::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double kRepuls, double lambda)
{
    double a = ScalarProductStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, b1->coordinate, b1->force, b2->coordinate, b2->force, lambda);
    double b = ScalarProductStretched(b3->coordinate, b3->force, b4->coordinate, b4->force, b3->coordinate, b3->force, b4->coordinate, b4->force, lambda);
    double c = ScalarProductStretched(b3->coordinate, b3->force, b1->coordinate, b1->force, b3->coordinate, b3->force, b1->coordinate, b1->force, lambda);
    double d = ScalarProductStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, b3->coordinate, b3->force, b4->coordinate, b4->force, lambda);
    double e = ScalarProductStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, b3->coordinate, b3->force, b1->coordinate, b1->force, lambda);
    double f = ScalarProductStretched(b3->coordinate, b3->force, b4->coordinate, b4->force, b3->coordinate, b3->force, b1->coordinate, b1->force, lambda);
    
    
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
    
    return 0.5*kRepuls/JJ*( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4 );
    
    
    
}

void CylinderExclVolRepulsion::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double kRepuls)
{
    double a = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b1->coordinateAux, b2->coordinateAux);
    double b = ScalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, b4->coordinateAux);
    double c = ScalarProduct(b3->coordinateAux, b1->coordinateAux, b3->coordinateAux, b1->coordinateAux);
    double d = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b3->coordinateAux, b4->coordinateAux);
    double e = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b3->coordinateAux, b1->coordinateAux);
    double f = ScalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, b1->coordinateAux);
    
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
    double A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) - (A2*CC*(a + e))/(AA*AA*AA);
    double A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
    double A14 = (A2*CC)/(AA*AA);
    
    double E11 = ATG2/EE;
    double E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) - (ATG2*GG)/(EE*EE);
    double E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
    double E14 = (E2*GG)/(EE*EE);
    
    double B11 = ATG3/BB;
    double B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*f)/(BB*BB*BB) + (B1*DD*(-b + f))/(BB*BB*BB);
    double B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
    double B14 = (B1*DD)/(BB*BB);
    
    double F11 = ATG4/FF;
    double F12 = (F2*(-d - f)*HH)/(FF*FF*FF) + (F1*(-b + d + f)*HH)/(FF*FF*FF) - (ATG4*HH)/(FF*FF);
    double F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
    double F14 = (F1*HH)/(FF*FF);

    //i1:(DIVIDING BY TWO)
    b1->force[0] +=  - 0.5*invJJ*( (b2->coordinate[0] - b1->coordinate[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinate[0] - b3->coordinate[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinate[0] - b3->coordinate[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    b1->force[1] +=  - 0.5*invJJ*( (b2->coordinate[1] - b1->coordinate[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinate[1] - b3->coordinate[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinate[1] - b3->coordinate[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    
    b1->force[2] +=  - 0.5*invJJ*( (b2->coordinate[2] - b1->coordinate[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinate[2] - b3->coordinate[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinate[2] - b3->coordinate[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    
    //i2:(DIVIDING BY TWO)
    b2->force[0] +=  - invJJ*( (b2->coordinate[0] - b1->coordinate[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinate[0] - b3->coordinate[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinate[0] - b3->coordinate[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->force[1] += - invJJ*( (b2->coordinate[1] - b1->coordinate[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinate[1] - b3->coordinate[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinate[1] - b3->coordinate[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->force[2] += - invJJ*( (b2->coordinate[2] - b1->coordinate[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinate[2] - b3->coordinate[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinate[2] - b3->coordinate[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    //j1:(DIVIDING BY TWO)
    b3->force[0] +=  - 0.5*invJJ*( (b2->coordinate[0] - b1->coordinate[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinate[0] - b3->coordinate[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinate[0] - b3->coordinate[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b3->force[1] +=  - 0.5*invJJ*( (b2->coordinate[1] - b1->coordinate[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinate[1] - b3->coordinate[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinate[1] - b3->coordinate[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
    
    b3->force[2] +=  - 0.5*invJJ*( (b2->coordinate[2] - b1->coordinate[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinate[2] - b3->coordinate[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinate[2] - b3->coordinate[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    //j2: (DIVIDING BY TWO)
    b4->force[0] +=  - invJJ*( 0.5*(b2->coordinate[0] - b1->coordinate[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinate[0] - b3->coordinate[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinate[0] - b3->coordinate[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) )  ;
    
    b4->force[1] +=  - invJJ*( 0.5*(b2->coordinate[1] - b1->coordinate[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinate[1] - b3->coordinate[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinate[1] - b3->coordinate[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->force[2] +=  - invJJ*( 0.5*(b2->coordinate[2] - b1->coordinate[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinate[2] - b3->coordinate[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinate[2] - b3->coordinate[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;

}

void CylinderExclVolRepulsion::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double kRepuls)
{
    double a = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b1->coordinateAux, b2->coordinateAux);
    double b = ScalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, b4->coordinateAux);
    double c = ScalarProduct(b3->coordinateAux, b1->coordinateAux, b3->coordinateAux, b1->coordinateAux);
    double d = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b3->coordinateAux, b4->coordinateAux);
    double e = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b3->coordinateAux, b1->coordinateAux);
    double f = ScalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, b1->coordinateAux);
    
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
    double A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) - (A2*CC*(a + e))/(AA*AA*AA);
    double A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
    double A14 = (A2*CC)/(AA*AA);
    
    double E11 = ATG2/EE;
    double E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) - (ATG2*GG)/(EE*EE);
    double E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
    double E14 = (E2*GG)/(EE*EE);
    
    double B11 = ATG3/BB;
    double B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*f)/(BB*BB*BB) + (B1*DD*(-b + f))/(BB*BB*BB);
    double B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
    double B14 = (B1*DD)/(BB*BB);
    
    double F11 = ATG4/FF;
    double F12 = (F2*(-d - f)*HH)/(FF*FF*FF) + (F1*(-b + d + f)*HH)/(FF*FF*FF) - (ATG4*HH)/(FF*FF);
    double F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
    double F14 = (F1*HH)/(FF*FF);
    
    //(DIVIDING BY TWO)
    b1->forceAux[0] +=  - 0.5*invJJ*( (b2->coordinateAux[0] - b1->coordinateAux[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinateAux[0] - b3->coordinateAux[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinateAux[0] - b3->coordinateAux[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
     //(DIVIDING BY TWO)
    b1->forceAux[1] +=  - 0.5*invJJ*( (b2->coordinateAux[1] - b1->coordinateAux[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinateAux[1] - b3->coordinateAux[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinateAux[1] - b3->coordinateAux[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
     //(DIVIDING BY TWO)
    b1->forceAux[2] +=  - 0.5*invJJ*( (b2->coordinateAux[2] - b1->coordinateAux[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (b4->coordinateAux[2] - b3->coordinateAux[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (b1->coordinateAux[2] - b3->coordinateAux[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
     //(DIVIDING BY TWO)
    b2->forceAux[0] +=  - invJJ*( (b2->coordinateAux[0] - b1->coordinateAux[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinateAux[0] - b3->coordinateAux[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinateAux[0] - b3->coordinateAux[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->forceAux[1] += - invJJ*( (b2->coordinateAux[1] - b1->coordinateAux[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinateAux[1] - b3->coordinateAux[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinateAux[1] - b3->coordinateAux[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
    b2->forceAux[2] += - invJJ*( (b2->coordinateAux[2] - b1->coordinateAux[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(b4->coordinateAux[2] - b3->coordinateAux[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(b1->coordinateAux[2] - b3->coordinateAux[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) );
    
     //(DIVIDING BY TWO)
    b3->forceAux[0] +=  - 0.5*invJJ*( (b2->coordinateAux[0] - b1->coordinateAux[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinateAux[0] - b3->coordinateAux[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinateAux[0] - b3->coordinateAux[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    b3->forceAux[1] +=  - 0.5*invJJ*( (b2->coordinateAux[1] - b1->coordinateAux[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinateAux[1] - b3->coordinateAux[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinateAux[1] - b3->coordinateAux[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) / 2;
    
    b3->forceAux[2] +=  - 0.5*invJJ*( (b2->coordinateAux[2] - b1->coordinateAux[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (b4->coordinateAux[2] - b3->coordinateAux[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (b1->coordinateAux[2] - b3->coordinateAux[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
     //(DIVIDING BY TWO)
    b4->forceAux[0] +=  - invJJ*( 0.5*(b2->coordinateAux[0] - b1->coordinateAux[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinateAux[0] - b3->coordinateAux[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinateAux[0] - b3->coordinateAux[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->forceAux[1] +=  - invJJ*( 0.5*(b2->coordinateAux[1] - b1->coordinateAux[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinateAux[1] - b3->coordinateAux[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinateAux[1] - b3->coordinateAux[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    b4->forceAux[2] +=  - invJJ*( 0.5*(b2->coordinateAux[2] - b1->coordinateAux[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (b4->coordinateAux[2] - b3->coordinateAux[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(b1->coordinateAux[2] - b3->coordinateAux[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) );
}

