
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

#include <cmath>
#include <vector>
#include <math.h>

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
            Closes distance is given by l = l0 + sc*l1 - tc*l2, where tc and sc are 
            numbers between[0,1], so points are withing the segments. Algorithm is the 
            following: find closes points between infinite lines and then move tc 
            and sc to 0 or 1 if they are out of [0,1].
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
    
    vector<double> movePointOutOfPlane(const vector<double>& p1,
                                       const vector<double>& p2,
                                       const vector<double>& p3,
                                       const vector<double>& p4,
                                       int i, double d) {
        vector<double> v;
        vector<double> v1;
        
        //plane
        v.push_back( (p2[2]-p1[2])*(p4[3]-p3[3]) -  (p2[3]-p1[3])*(p4[2]-p3[2]) );
        v.push_back( (p2[3]-p1[3])*(p4[1]-p3[1]) -  (p2[1]-p1[1])*(p4[3]-p3[3]) );
        v.push_back( (p2[1]-p1[1])*(p4[2]-p3[2]) -  (p2[2]-p1[2])*(p4[1]-p3[1]) );
        
        double norm = sqrt( v[1]*v[1] + v[2]*v[2] + v[3]*v[3] );
        
        v1.push_back(v[0]/norm);
        v1.push_back(v[1]/norm);
        v1.push_back(v[2]/norm);
        
        //move bead 1
        if (i == 1){
            vector<double> newP1;
            newP1.push_back(p1[0] + v1[0]*d);
            newP1.push_back(p1[1] + v1[1]*d);
            newP1.push_back(p1[2] + v1[2]*d);
            return newP1;
        }
        
        //move bead 2
        else if (i == 2){
            vector<double> newP2;
            newP2.push_back(p2[0] + v1[0]*d);
            newP2.push_back(p2[1] + v1[1]*d);
            newP2.push_back(p2[2] + v1[2]*d);
            return newP2;
        }
        
        //move bead 3
        else if (i == 3){
            vector<double> newP3;
            newP3.push_back(p3[0] + v1[0]*d);
            newP3.push_back(p3[1] + v1[1]*d);
            newP3.push_back(p3[2] + v1[2]*d);
            return newP3;
        }
        
        //move bead 4
        else {
            vector<double> newP4;
            newP4.push_back(p4[0] + v1[0]*d);
            newP4.push_back(p4[1] + v1[1]*d);
            newP4.push_back(p4[2] + v1[2]*d);
            return newP4;
        }
    }
    
    
    tuple<vector<double>, vector<double>> branchProjection(const vector<double>& n,
                                                           const vector<double>& p,
                                                           double l, double m, double theta){
        //get random permutation from p
        vector<double> r = {p[0] + randomDouble(-1, 1),
                            p[1] + randomDouble(-1, 1),
                            p[2] + randomDouble(-1, 1)};
        
        //construct vector z which is r-p
        auto z = twoPointDirection(p, r);
        
        //construct u and v, which creates an orthogonal set n, u, v
        auto u = crossProduct(n, z);
        auto v = crossProduct(n, u);
        
        normalize(u); normalize(v);
        
        //find random point on circle defining the branching point
        double thetaRandom = randomDouble(0, 2*M_PI);
        vector<double> bp1;
        bp1.push_back(p[0] + l * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp1.push_back(p[1] + l * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp1.push_back(p[2] + l * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));
        
        //now find the second point
        vector<double> newP;
        double dist = m * cos(theta);
        newP.push_back(p[0] + n[0] * dist);
        newP.push_back(p[1] + n[1] * dist);
        newP.push_back(p[2] + n[2] * dist);
        double newL = (l + m * sin(theta));
        
        vector<double> bp2;
        bp2.push_back(newP[0] + newL * ( u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp2.push_back(newP[1] + newL * ( u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp2.push_back(newP[2] + newL * ( u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));
            
        //get direction
        auto direction = twoPointDirection(bp1, bp2);
        
        return tuple<vector<double>, vector<double>>(direction, bp1);
    }
}




