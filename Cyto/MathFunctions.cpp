//
//  MathFunctions.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MathFunctions.h"
#include <cmath>
#include <vector>

using namespace std;

namespace mathfunc {
    
    
    /*! This is a collection of diferent type of mathematical functions which are use through entire program;
     */
    
    /// A function which calculates a nuber — distance between two points seted up by two vectors v1 = {x1, y1, z1} and v1 = {x2, y2, z2};
    double TwoPointDistance(const std::vector<double>& v1, const std::vector<double>& v2) {
        
        double l= sqrt( (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
        
        //    cout<<"Dist = "<<l<<endl;
        return l;
        
    }
    
    /// A  auxilury function which calculates a nuber — distance between two points seted up by two vectors v1-d*p1 = {x1-d*px1, y1-d*py1, z1-d*pz1} and v2-d*p2 = {x2-d*px2, y2-d*py2, z2-d*pz2}; needed for Golden Section minimization;
    
    double TwoPointDistanceStretched(const std::vector<double>& v1, const std::vector<double>& p1, const std::vector<double>& v2, const std::vector<double>& p2, double d ){
        
        double l= sqrt( ((v2[0] + d*p2[0])-(v1[0] + d*p1[0]))*((v2[0] + d*p2[0])-(v1[0] + d*p1[0])) + ((v2[1] + d*p2[1])-(v1[1] + d*p1[1]))*((v2[1] + d*p2[1])-(v1[1] + d*p1[1])) + ((v2[2] + d*p2[2])-(v1[2] + d*p1[2]))*((v2[2] + d*p2[2])-(v1[2] + d*p1[2])));
        
        //    cout<<"Dist_str = "<<l<<endl;
        return l;
    }
    
    
    vector<double> TwoPointDirection(const std::vector<double>& v1, const std::vector<double>& v2) {
        vector<double> tau (3, 0);
        double invD = 1/TwoPointDistance(v1, v2);
        tau[0] = invD * ( v2[0] - v1[0] );
        tau[1] = invD * ( v2[1] - v1[1] );
        tau[2] = invD * ( v2[2] - v1[2] );
        return tau;
        
    }
    ///A function which returns coordinate for a next point projection based on initial coordinates and a diraction;
    vector<double> NextPointProjection(const std::vector<double>& coordinate, double d, const std::vector<double>& tau){
        
        vector<double> v;
        v.push_back(coordinate[0] + d * tau[0]);
        v.push_back(coordinate[1] + d * tau[1]);
        v.push_back(coordinate[2] + d * tau[2]);
        
        return v;
        
        
    }
    
    double ScalarProduct(const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& v3, const std::vector<double>& v4){ return ( (v2[0]-v1[0])*(v4[0]-v3[0]) + (v2[1]-v1[1])*(v4[1]-v3[1]) +(v2[2]-v1[2])*(v4[2]-v3[2]) );}
    
    double ScalarProductStretched(const std::vector<double>& v1, const std::vector<double>& p1, const std::vector<double>& v2, const std::vector<double>& p2, const std::vector<double>& v3, const std::vector<double>& p3, const std::vector<double>& v4, const std::vector<double>& p4, double d){
        
        double xx = ( (v2[0] - d*p2[0])-(v1[0] - d*p1[0]) )*( (v4[0] - d*p4[0])-(v3[0] - d*p3[0]) );
        double yy = ( (v2[1] - d*p2[1])-(v1[1] - d*p1[1]) )*( (v4[1] - d*p4[1])-(v3[1] - d*p3[1]) );
        double zz = ( (v2[2] - d*p2[2])-(v1[2] - d*p1[2]) )*( (v4[2] - d*p4[2])-(v3[2] - d*p3[2]) );
        
        return xx + yy + zz;
        
    }
    
    std::vector<double> MidPointCoordinate(const std::vector<double>& v1, const std::vector<double>& v2, double alpha){
        std::vector<double> v;
        
        v.push_back( v1[0]*(1.0 - alpha) + alpha*v2[0] );
        v.push_back( v1[1]*(1.0 - alpha) + alpha*v2[1] );
        v.push_back( v1[2]*(1.0 - alpha) + alpha*v2[2] );
        
        return v;
    }
    
    
    std::vector<double> MidPointCoordinateStretched(const std::vector<double>& v1, const std::vector<double>& p1, const std::vector<double>& v2, const std::vector<double>& p2, double alpha , double d){
        
        std::vector<double> v;
        
        v.push_back( (v1[0] - d*p1[0])*(1.0 - alpha) + alpha*(v2[0] - d*p2[0]) );
        v.push_back( (v1[1] - d*p1[1])*(1.0 - alpha) + alpha*(v2[1] - d*p2[1]) );
        v.push_back( (v1[2] - d*p1[2])*(1.0 - alpha) + alpha*(v2[2] - d*p2[2]) );
        return v;
    }
    
    double TwoSegmentDistance(const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& v3, const std::vector<double>& v4){
        
        
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
            
            lines are given by a parametric equation v1 + s*(v2-v1) and v3 + t*(v4-v3). Closes distance is given by l = l0 + sc*l1 - tc*l2, where tc and sc are numbers between[0,1], so points are withing the segments. Algorythm is the following: find closes points between infinite lines and then move tc and sc to 0 or 1 if they are out of [0,1].
 
*/
        
        double SMALL_NUM = 0.0001;
        double a = ScalarProduct(v1, v2, v1, v2);
        double b = ScalarProduct(v1, v2, v3, v4);
        double c = ScalarProduct(v3, v4, v3, v4);
        double d = ScalarProduct(v1, v2, v1, v3);
        double e = ScalarProduct(v3, v4, v1, v3);
        double f = ScalarProduct(v1, v3, v1, v3);
        
        
        
        double D = a*c - b*b;
        double sc, sN, tc, tN;
        double sD = D;
        double tD = D;
        
        
        
        // compute the line parameters of the two closest points
        if (D < SMALL_NUM) { // the lines are almost parallel
            sN = 0.0;         // force using point v1 on segment l1
            sD = 1.0;         // to prevent possible division by 0.0 later
            tN = e;
            tD = c;
            }
        else {                 // get the closest points on the infinite lines
            sN = (b*e - c*d);
            tN = (a*e - b*d);
            if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
                sN = 0.0;
                tN = e;
                tD = c;
            }
            else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }
            
        if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
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
        else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
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
        
        sc = ( abs( sN ) < SMALL_NUM ? 0.0 : sN / sD);
        tc = ( abs( tN ) < SMALL_NUM ? 0.0 : tN / tD);
    
       
        
        // Calculate |(l0+sc*l1-tc*l2)|:
        return sqrt( f + sc*sc*a +tc*tc*c + 2*(sc*d - tc*e - tc*sc*b) );
        
  
    }
    
}
