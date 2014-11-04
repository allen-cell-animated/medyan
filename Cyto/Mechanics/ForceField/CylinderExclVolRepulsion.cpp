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
#include <math.h>


using namespace std;
using namespace mathfunc;


double CylinderExclVolRepulsion::Energy(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double kRepuls)
{
    double a = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate);
    double b = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb4->coordinate);
    double c = ScalarProduct(pb3->coordinate, pb1->coordinate, pb3->coordinate, pb1->coordinate);
    double d = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb4->coordinate);
    double e = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb1->coordinate);
    double f = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb1->coordinate);
    
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

double CylinderExclVolRepulsion::Energy(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double kRepuls, double lambda)
{
    double a = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, lambda);
    double b = ScalarProductStretched(pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, lambda);
    double c = ScalarProductStretched(pb3->coordinate, pb3->force, pb1->coordinate, pb1->force, pb3->coordinate, pb3->force, pb1->coordinate, pb1->force, lambda);
    double d = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, lambda);
    double e = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb3->coordinate, pb3->force, pb1->coordinate, pb1->force, lambda);
    double f = ScalarProductStretched(pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, pb3->coordinate, pb3->force, pb1->coordinate, pb1->force, lambda);
    
    
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

void CylinderExclVolRepulsion::Forces(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double kRepuls)
{
    double a = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate);
    double b = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb4->coordinate);
    double c = ScalarProduct(pb3->coordinate, pb1->coordinate, pb3->coordinate, pb1->coordinate);
    double d = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb4->coordinate);
    double e = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb1->coordinate);
    double f = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb1->coordinate);
    
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
    double A12 = -((ATG1*CC)/AA*AA) + (A1*CC*e)/AA*AA*AA - (A2*CC*(a + e))/AA*AA*AA;
    double A13 = -((A1*CC)/AA*AA) + (A2*CC)/AA*AA;
    double A14 = (A2*CC)/AA;
    
    double E11 = ATG2/EE;
    double E12 = (E2*(-a + d - e)*GG)/EE*EE*EE + (E1*(-d + e)*GG)/EE*EE*EE - (ATG2*GG)/EE*EE;
    double E13 = -((E1*GG)/EE*EE) + (E2*GG)/EE*EE;
    double E14 = (E2*GG)/EE*EE;
    
    double B11 = ATG3/BB;
    double B12 = -((ATG3*DD)/BB*BB) - (B2*DD*f)/BB*BB*BB + (B1*DD*(-b + f))/BB*BB*BB;
    double B13 = -((B1*DD)/BB*BB) + (B2*DD)/BB*BB;
    double B14 = (B1*DD)/BB*BB;
    
    double F11 = ATG4/FF;
    double F12 = (F2*(-d - f)*HH)/FF*FF*FF + (F1*(-b + d + f)*HH)/FF*FF - (ATG4*HH)/FF*FF;
    double F13 = -((F1*HH)/FF*FF) + (F2*HH)/FF*FF;
    double F14 = (F1*HH)/FF*FF;
    
    
    pb1->force[0] +=  - 0.5*invJJ*( (pb2->coordinate[0] - pb1->coordinate[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f -
        2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (pb4->coordinate[0] - pb3->coordinate[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (pb1->coordinate[0] - pb3->coordinate[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
    
   
    pb1->force[1] +=  - 0.5*invJJ*( (pb2->coordinate[1] - pb1->coordinate[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (pb4->coordinate[1] - pb3->coordinate[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (pb1->coordinate[1] - pb3->coordinate[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );


    pb1->force[2] +=  - 0.5*invJJ*( (pb2->coordinate[2] - pb1->coordinate[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (pb4->coordinate[2] - pb3->coordinate[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (pb1->coordinate[2] - pb3->coordinate[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

    

    pb2->force[0] +=  - invJJ*( (pb2->coordinate[0] - pb1->coordinate[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(pb4->coordinate[0] - pb3->coordinate[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(pb1->coordinate[0] - pb3->coordinate[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) ) ;
    
    pb2->force[1] += - invJJ*( (pb2->coordinate[1] - pb1->coordinate[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(pb4->coordinate[1] - pb3->coordinate[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(pb1->coordinate[1] - pb3->coordinate[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) ) ;
    
    pb2->force[2] += - invJJ*( (pb2->coordinate[2] - pb1->coordinate[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(pb4->coordinate[2] - pb3->coordinate[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(pb1->coordinate[2] - pb3->coordinate[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) ) ;
    
    
    pb3->force[0] +=  - 0.5*invJJ*( (pb2->coordinate[0] - pb1->coordinate[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (pb4->coordinate[0] - pb3->coordinate[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (pb1->coordinate[0] - pb3->coordinate[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );

    pb3->force[1] +=  - 0.5*invJJ*( (pb2->coordinate[1] - pb1->coordinate[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (pb4->coordinate[1] - pb3->coordinate[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (pb1->coordinate[1] - pb3->coordinate[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );

    pb3->force[2] +=  - 0.5*invJJ*( (pb2->coordinate[2] - pb1->coordinate[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (pb4->coordinate[2] - pb3->coordinate[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (pb1->coordinate[2] - pb3->coordinate[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );


    
    pb4->force[0] +=  - invJJ*( 0.5*(pb2->coordinate[0] - pb1->coordinate[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (pb4->coordinate[0] - pb3->coordinate[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(pb1->coordinate[0] - pb3->coordinate[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;



}

void CylinderExclVolRepulsion::ForcesAux(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double kRepuls)
{
    double a = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate);
    double b = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb4->coordinate);
    double c = ScalarProduct(pb3->coordinate, pb1->coordinate, pb3->coordinate, pb1->coordinate);
    double d = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb4->coordinate);
    double e = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb1->coordinate);
    double f = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb1->coordinate);
    
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
    double A12 = -((ATG1*CC)/AA*AA) + (A1*CC*e)/AA*AA*AA - (A2*CC*(a + e))/AA*AA*AA;
    double A13 = -((A1*CC)/AA*AA) + (A2*CC)/AA*AA;
    double A14 = (A2*CC)/AA;
    
    double E11 = ATG2/EE;
    double E12 = (E2*(-a + d - e)*GG)/EE*EE*EE + (E1*(-d + e)*GG)/EE*EE*EE - (ATG2*GG)/EE*EE;
    double E13 = -((E1*GG)/EE*EE) + (E2*GG)/EE*EE;
    double E14 = (E2*GG)/EE*EE;
    
    double B11 = ATG3/BB;
    double B12 = -((ATG3*DD)/BB*BB) - (B2*DD*f)/BB*BB*BB + (B1*DD*(-b + f))/BB*BB*BB;
    double B13 = -((B1*DD)/BB*BB) + (B2*DD)/BB*BB;
    double B14 = (B1*DD)/BB*BB;
    
    double F11 = ATG4/FF;
    double F12 = (F2*(-d - f)*HH)/FF*FF*FF + (F1*(-b + d + f)*HH)/FF*FF - (ATG4*HH)/FF*FF;
    double F13 = -((F1*HH)/FF*FF) + (F2*HH)/FF*FF;
    double F14 = (F1*HH)/FF*FF;
    
    
    pb1->forceAux[0] +=  - 0.5*invJJ*( (pb2->coordinate[0] - pb1->coordinate[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (pb4->coordinate[0] - pb3->coordinate[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (pb1->coordinate[0] - pb3->coordinate[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
    
    
    pb1->forceAux[1] +=  - 0.5*invJJ*( (pb2->coordinate[1] - pb1->coordinate[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (pb4->coordinate[1] - pb3->coordinate[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (pb1->coordinate[1] - pb3->coordinate[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    pb1->forceAux[2] +=  - 0.5*invJJ*( (pb2->coordinate[2] - pb1->coordinate[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*f))/(2*EE) - A11*f + E11*f - 2*U*f*f + (F12*b)/(2*FF)) ) + (pb4->coordinate[2] - pb3->coordinate[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*f - F11*f - 2*U*a*f - 4*U*e*f + 2*U*(d*e - a*f) - (B12*f)/BB) +  (pb1->coordinate[2] - pb3->coordinate[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*f + 2*U*(b*e - d*f) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    
    pb2->forceAux[0] +=  - invJJ*( (pb2->coordinate[0] - pb1->coordinate[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(pb4->coordinate[0] - pb3->coordinate[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(pb1->coordinate[0] - pb3->coordinate[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) ) ;
    
    pb2->forceAux[1] += - invJJ*( (pb2->coordinate[1] - pb1->coordinate[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(pb4->coordinate[1] - pb3->coordinate[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(pb1->coordinate[1] - pb3->coordinate[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) ) ;
    
    pb2->forceAux[2] += - invJJ*( (pb2->coordinate[2] - pb1->coordinate[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*f))/(2*EE)-A11*f+E11*f-2*U*f*f+(F12*b)/(2*FF) ) + 0.5*(pb4->coordinate[2] - pb3->coordinate[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f + 4*U*e*f - (F12*(d + f))/FF)  + 0.5*(pb1->coordinate[2] - pb3->coordinate[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*f - 2*U*(b*e - d*f) + (F12*b)/FF) ) ;
    
    
    pb3->forceAux[0] +=  - 0.5*invJJ*( (pb2->coordinate[0] - pb1->coordinate[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (pb4->coordinate[0] - pb3->coordinate[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (pb1->coordinate[0] - pb3->coordinate[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    pb3->forceAux[1] +=  - 0.5*invJJ*( (pb2->coordinate[1] - pb1->coordinate[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (pb4->coordinate[1] - pb3->coordinate[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (pb1->coordinate[1] - pb3->coordinate[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    pb3->forceAux[2] +=  - 0.5*invJJ*( (pb2->coordinate[2] - pb1->coordinate[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*f - F11*f - 2*U*d*f - 4*U*e*f + 2*U*(b*e - d*f) - (F12*b)/FF + (F12*(d + f))/FF) + (pb4->coordinate[2] - pb3->coordinate[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (pb1->coordinate[2] - pb3->coordinate[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*f - 2*U*(d*e - a*f) + (B12*f)/BB + (F12*(d + f))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
    
    
    
    pb4->forceAux[0] +=  - invJJ*( 0.5*(pb2->coordinate[0] - pb1->coordinate[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*f + F11*f +4*U*e*f - (F12*(d + f))/FF ) + (pb4->coordinate[0] - pb3->coordinate[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(pb1->coordinate[0] - pb3->coordinate[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*f + 2*U*(d*e - a*f) - (B12*f)/BB - (F12*(d + f))/FF) ) ;
    
    
    
}

