//
//  VolRepulsionForce.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 7/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "VolumInter.h"
#include "Bead.h"
#include <cmath>

using namespace std;
using namespace mathfunc;

//void VolumInter::ForceRepulsion(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4 ){
//     double a, b, c, d, e, f, A, B, C, D, E, F, G, H, J;
    
     //Segment one: between bead 2 and 1, coordinates x1-----x2; Segment two: between bead 4 and bead 3, coordinates: y1------y2.
     //x = x2 - x1, y = y2 - y1, z = y1 - x1;
     
//     a = (x,x)
 //    b = (y,y)
//     c = (z,z)
//     d = (x,y)
//     e = (x,z)
//     f = (y,z)
     
    
//    a = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate);
//    b = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb4->coordinate);
//    c = ScalarProduct(pb1->coordinate, pb3->coordinate, pb1->coordinate, pb3->coordinate);
    
//    d = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb4->coordinate);
//    e = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb3->coordinate);
//    f = ScalarProduct(pb3->coordinate, pb4->coordinate, pb1->coordinate, pb3->coordinate);
    
//    A = sqrt(a*c - e*e);
//    B = sqrt(b*c - f*f);
    
//    C = d*e - a*f;
//    D = b*e - d*f;
    
//    E = sqrt( a*(b+c-2*f) - (d-e)*(d-e) );
//    F = sqrt( b*(a+c-2*e) - (d+f)*(d+f) );
    
//    G = d*d  - a*b - C;
//    H = G + C - D;
//    J = c*(G + C) + e*D - f*C;
    
//    double ATG1 = atan((a+e)/A)  - atan(e/A);
//    double ATG2 = atan((a+e-d)/E)  - atan((e-d)/E);
//    double ATG3 = atan((f)/B)  - atan((f-b)/B);
//    double ATG4 = atan((d+f)/F)  - atan((d+f-b)/F);
    
//    double A1 = A*A/(A*A + e*e);
//    double A2 = A*A/(A*A + (a+e)*(a+e));
    
//    double E1 = E*E/(E*E + (a+e-d)*(a+e-d));
//    double E2 = E*E/(E*E + (e-d)*(e-d));
    
//    double B1 = B*B/(B*B + f*f);
//    double B2 = B*B/(B*B + (f-b)*(f-b));
    
//    double F1 = F*F/(F*F + (d+f)*(d+f));
//    double F2 = B*B/(B*B + (d+f-b)*(d+f-b));
    
    
    
//}