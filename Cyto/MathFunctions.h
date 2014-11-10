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

#include "common.h"

namespace mathfunc {
    
    
    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3);
    double TwoPointDistance(const vector<double>& v1, const vector<double>& v2);
    
    /// Compute distance between two points with coordinates: (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z). Needed for Golden Section Minimization;
    
    double TwoPointDistanceStretched(const vector<double>& v1, const vector<double>& p1,
                                     const vector<double>& v2, const vector<double>& p2, double d );
    
    ///Calculates a normal to a line starting at (x1,y1,z1) and ending at (x2,y2,z2);
    vector<double> TwoPointDirection(const vector<double>& v1,
                                          const vector<double>& v2);
    
    ///Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and (x4-x3,y4-y3,z4-z3);
    double ScalarProduct(const vector<double>& v1, const vector<double>& v2,
                         const vector<double>& v3, const vector<double>& v4);
    
    ///Same as the above but with x-d*p coordinates;
    double ScalarProductStretched(const vector<double>& v1, const vector<double>& p1,
                                  const vector<double>& v2, const vector<double>& p2,
                                  const vector<double>& v3, const vector<double>& p3,
                                  const vector<double>& v4, const vector<double>& p4, double d);
    
    vector<double> NextPointProjection(const vector<double>&, double, const vector<double>&);
    
    /// Motor related function, returns coordinates of a point v located on a line between v1 and v2. |v-v1|/|v2-v| = alpha.
    vector<double> MidPointCoordinate(const vector<double>&, const vector<double>&, double);
    
    vector<double> MidPointCoordinateStretched(const vector<double>& v1, const vector<double>& p1,
                                               const vector<double>& v2, const vector<double>& p2,
                                               double alpha , double d);
    
    /// Function to calculate a diatance between two segment. Used in excluded volume interaction calculations.
    double TwoSegmentDistance(const vector<double>& v1, const vector<double>& v2,
                              const vector<double>& v3, const vector<double>& v4);
    
}



#endif /* defined(__CytoMech__MathFunctions__) */
