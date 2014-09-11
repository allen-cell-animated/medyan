//
//  MathFunctions.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MathFunctions__
#define __CytoMech__MathFunctions__

#include <iostream>
#include <cmath>
#include <vector>

namespace mathfunc {
    
    
    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3);
    double TwoPointDistance(std::vector<double> v1, std::vector<double> v2);
    
    /// Compute distance between two points with coordinates: (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z). Needed for Golden Section Minimization;
    
    double TwoPointDistanceStretched(std::vector<double> v1, std::vector<double> p1, std::vector<double> v2, std::vector<double> p2, double d );
    
    ///Calculates a normal to a line starting at (x1,y1,z1) and ending at (x2,y2,z3);
    std::vector<double> TwoPointDirection(std::vector<double> v1, std::vector<double> v2);
    
    ///Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and (x4-x3,y4-y3,z4-z3);
    double ScalarProduct(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> v4);
    
    ///Same as the above but with x-d*p coordinates;
    double ScalarProductStretched(std::vector<double> v1, std::vector<double> p1, std::vector<double> v2, std::vector<double> p2, std::vector<double> v3, std::vector<double> p3, std::vector<double> v4, std::vector<double> p4, double d);
    
    std::vector<double> NextPointProjection(std::vector<double>, double, std::vector<double>);
    
    /// Motor related function, returns coordinates of a point v located on a line between v1 and v2. |v-v1|/|v2-v| = alpha.
    std::vector<double> MidPointCoordinate(std::vector<double>, std::vector<double>, double);
    std::vector<double> MidPointCoordinateStretched(std::vector<double> v1, std::vector<double> p1, std::vector<double> v2, std::vector<double> p2, double alpha , double d);
    
    /// Function to calculate a diatance between two segment. Used in excluded volume interaction calculations.
    double TwoSegmentDistance(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> v4);
    
}



#endif /* defined(__CytoMech__MathFunctions__) */
