
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




