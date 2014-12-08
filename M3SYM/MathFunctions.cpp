
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include <cmath>
#include <vector>

#include "MathFunctions.h"

namespace mathfunc {
    
    double twoSegmentDistance(const vector<double>& v1, const vector<double>& v2,
                              const vector<double>& v3, const vector<double>& v4){
        
/*!         
            segment 1:  l1 = v2 - v1;
            segment 2:  l2 = v3 - v4;
            distance betwin begining of l1 nd l2:   l0 =  v3 - v1;
            a = (l1,l1);
            b = (l1,l2);
            c = (ls,vls);
            d = (l1,l0);
            e = (l2,l0);
            f = (l0,l0);
            D = a*c - b*b;
            
            lines are given by a parametric equation v1 + s*(v2-v1) and v3 + t*(v4-v3). 
            Closes distance is given by l = l0 + sc*l1 - tc*l2, where tc and sc are numbers between[0,1], 
            so points are withing the segments. Algorithm is the following: find closes points between infinite 
            lines and then move tc and sc to 0 or 1 if they are out of [0,1].
*/
        
        double SMALL_NUM = 0.0001;
        double a = scalarProduct(v1, v2, v1, v2);
        double b = scalarProduct(v1, v2, v3, v4);
        double c = scalarProduct(v3, v4, v3, v4);
        double d = scalarProduct(v1, v2, v1, v3);
        double e = scalarProduct(v3, v4, v1, v3);
        double f = scalarProduct(v1, v3, v1, v3);
        
        double D = a*c - b*b;
        double sc, sN, tc, tN;
        double sD = D;
        double tD = D;
        
        // compute the line parameters of the two closest points
        
        // the lines are almost parallel
        if (D < SMALL_NUM) {
            sN = 0.0;         // force using point v1 on segment l1
            sD = 1.0;         // to prevent possible division by 0.0 later
            tN = e;
            tD = c;
        }
        // get the closest points on the infinite lines
        else {
            sN = (b*e - c*d);
            tN = (a*e - b*d);
            
            // sc < 0 => the s=0 edge is visible
            if (sN < 0.0) {
                sN = 0.0;
                tN = e;
                tD = c;
            }
            // sc > 1  => the s=1 edge is visible
            else if (sN > sD) {
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }
        // tc < 0 => the t=0 edge is visible
        if (tN < 0.0) {
            tN = 0.0;
            // recompute sc for this edge
            if (-d < 0.0)
                sN = 0.0;
            else if (-d > a)
                sN = sD;
            else {
                sN = -d;
                sD = a;
            }
        }
        // tc > 1  => the t=1 edge is visible
        else if (tN > tD) {
                tN = tD;
            // recompute sc for this edge
            if ((-d + b) < 0.0)
                sN = 0;
            else if ((-d + b) > a)
                sN = sD;
            else {
                sN = (-d +  b);
                sD = a;
            }
        }
        
        // Compute sc and tc
        sc = ( fabs( sN ) < SMALL_NUM ? 0.0 : sN / sD);
        tc = ( fabs( tN ) < SMALL_NUM ? 0.0 : tN / tD);
    
        // Calculate |(l0+sc*l1-tc*l2)|:
        return sqrt(f + sc*sc*a +tc*tc*c + 2*(sc*d - tc*e - tc*sc*b));
    }
}
