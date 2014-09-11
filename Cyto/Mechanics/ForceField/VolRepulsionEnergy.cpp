//
//  VolRepulsion.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 6/26/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "VolumInter.h"
#include "Bead.h"
#include <cmath>

using namespace std;
using namespace mathfunc;
//
//double VolumInter::EnergyRepulsion(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4)
//{
//    
//    double U = 0.0;
//    double a, b, c, d, e, f, A, B, C, D, E, F, G, H, J;
//    /*
//     Segment one: between bead 2 and 1, coordinates x1-----x2; Segment two: between bead 4 and bead 3, coordinates: y1------y2.
//     x = x2 - x1, y = y2 - y1, z = y1 - x1;
//     
//     a = (x,x)
//     b = (y,y)
//     c = (z,z)
//     d = (x,y)
//     e = (x,z)
//     f = (y,z)
//     */
//    
//    a = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate);
//    b = ScalarProduct(pb3->coordinate, pb4->coordinate, pb3->coordinate, pb4->coordinate);
//    c = ScalarProduct(pb1->coordinate, pb3->coordinate, pb1->coordinate, pb3->coordinate);
//    
//    d = ScalarProduct(pb1->coordinate, pb2->coordinate, pb3->coordinate, pb4->coordinate);
//    e = ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb3->coordinate);
//    f = ScalarProduct(pb3->coordinate, pb4->coordinate, pb1->coordinate, pb3->coordinate);
//    
//    A = sqrt(a*c - e*e);
//    B = sqrt(b*c - f*f);
//    
//    C = d*e - a*f;
//    D = b*e - d*f;
//    
//    E = sqrt( a*(b+c-2*f) - (d-e)*(d-e) );
//    F = sqrt( b*(a+c-2*e) - (d+f)*(d+f) );
//    
//    G = d*d  - a*b - C;
//    H = G + C - D;
//    J = c*(G + C) + e*D - f*C;
//    
//    U = _kRepuls* 0.5/J * ( C/A*( atan((a+e)/A)  - atan(e/A) ) + G/E*( atan((a+e-d)/E)  - atan((e-d)/E) ) + D/B*( atan((f)/B)  - atan((f-b)/B) ) + H/F*( atan((d+f)/F)  - atan((d+f-b)/F) ));
//    
//    return U;
//}
//
//// Compute energy with x -> x-p*f argement. For CG method.
//double VolumInter::EnergyRepulsion(Bead *pb1, Bead *pb2, Bead *pb3, Bead *pb4, double l)
//{
//    
//    double U = 0.0;
//    double a, b, c, d, e, f, A, B, C, D, E, F, G, H, J;
//    /*
//     Segment one: between bead 2 and 1, coordinates x1-----x2; Segment two: between bead 4 and bead 3, coordinates: y1------y2.
//     x = x2 - x1, y = y2 - y1, z = y1 - x1;
//     
//     a = (x,x)
//     b = (y,y)
//     c = (z,z)
//     d = (x,y)
//     e = (x,z)
//     f = (y,z)
//     */
//    
//    a = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, l);
//    b = ScalarProductStretched(pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, l);
//    c = ScalarProductStretched(pb1->coordinate, pb1->force, pb3->coordinate, pb3->force, pb1->coordinate, pb1->force, pb3->coordinate, pb3->force, l);
//    
//    d = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, l);
//    e = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb1->coordinate, pb1->force, pb3->coordinate, pb3->force, l);
//    f = ScalarProductStretched(pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, pb1->coordinate, pb1->force, pb3->coordinate, pb3->force, l);
//    
//    A = sqrt(a*c - e*e);
//    B = sqrt(b*c - f*f);
//    
//    C = d*e - a*f;
//    D = b*e - d*f;
//    
//    E = sqrt( a*(b+c-2*f) - (d-e)*(d-e) );
//    F = sqrt( b*(a+c-2*e) - (d+f)*(d+f) );
//    
//    G = d*d  - a*b - C;
//    H = G + C - D;
//    J = c*(G + C) + e*D - f*C;
//    
//    U = _kRepuls* 0.5/J * ( C/A*( atan((a+e)/A)  - atan(e/A) ) + G/E*( atan((a+e-d)/E)  - atan((e-d)/E) ) + D/B*( atan((f)/B)  - atan((f-b)/B) ) + H/F*( atan((d+f)/F)  - atan((d+f-b)/F) ));
//    
//    return U;
//}
//
