
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

#ifndef M3SYM_MathFunctions_h
#define M3SYM_MathFunctions_h

#include <cmath>
#include <vector>

#include "common.h"

/// @namespace mathfunc is used for the mathematics module for the entire codebase
/// mathfunc includes functions to calculate distances, products, and midpoints
namespace mathfunc {
    
    /// Normalize a vector
    inline void normalize(vector<double>& v) {
        
        double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }
    
    
    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3)
    inline double twoPointDistance(const vector<double>& v1, const vector<double>& v2) {
        return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) +
                    (v2[1]-v1[1])*(v2[1]-v1[1]) +
                    (v2[2]-v1[2])*(v2[2]-v1[2]));
    }
    
    /// Compute distance between two points with coordinates
    /// (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z)
    inline double twoPointDistanceStretched(const vector<double>& v1,
                                            const vector<double>& p1,
                                            const vector<double>& v2,
                                            const vector<double>& p2, double d){
        
        return sqrt(((v2[0] + d*p2[0])-(v1[0] + d*p1[0])) *
                    ((v2[0] + d*p2[0])-(v1[0] + d*p1[0])) +
                    ((v2[1] + d*p2[1])-(v1[1] + d*p1[1])) *
                    ((v2[1] + d*p2[1])-(v1[1] + d*p1[1])) +
                    ((v2[2] + d*p2[2])-(v1[2] + d*p1[2])) *
                    ((v2[2] + d*p2[2])-(v1[2] + d*p1[2])));
    }

    
    /// Calculates a normal to a line starting at (x1,y1,z1) and ending at (x2,y2,z2)
    inline vector<double> twoPointDirection(const vector<double>& v1,
                                            const vector<double>& v2) {
        vector<double> tau (3, 0);
        double invD = 1/twoPointDistance(v1, v2);
        tau[0] = invD * ( v2[0] - v1[0] );
        tau[1] = invD * ( v2[1] - v1[1] );
        tau[2] = invD * ( v2[2] - v1[2] );
        return tau;
    }
    
    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    inline double scalarProduct(const vector<double>& v1, const vector<double>& v2,
                                const vector<double>& v3, const vector<double>& v4) {
        return ((v2[0] - v1[0])*(v4[0] - v3[0]) +
                (v2[1] - v1[1])*(v4[1] - v3[1]) +
                (v2[2] - v1[2])*(v4[2] - v3[2]));
    }
    
    /// Scalar product of two vectors v1(x,y,z) and v2(x,y,z)
    inline double dotProduct(const vector<double>& v1, const vector<double>& v2) {
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }

    
    
    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3) but with x+d*p coordinates
    inline double scalarProductStretched(const vector<double>& v1,
                                         const vector<double>& p1,
                                         const vector<double>& v2,
                                         const vector<double>& p2,
                                         const vector<double>& v3,
                                         const vector<double>& p3,
                                         const vector<double>& v4,
                                         const vector<double>& p4, double d){
        
        double xx = ((v2[0] + d*p2[0]) - (v1[0] + d*p1[0]))*
                    ((v4[0] + d*p4[0]) - (v3[0] + d*p3[0]));
        double yy = ((v2[1] + d*p2[1]) - (v1[1] + d*p1[1]))*
                    ((v4[1] + d*p4[1]) - (v3[1] + d*p3[1]));
        double zz = ((v2[2] + d*p2[2]) - (v1[2] + d*p1[2]))*
                    ((v4[2] + d*p4[2]) - (v3[2] + d*p3[2]));
        return xx + yy + zz;
        
    }
    
    /// Scalar product of two vectors with coordinates: v1[z,y,z] + d*p1[x,y,z] and
    /// v2[x,y,z] + d*p2[x,y,z]
    inline double dotProductStretched(const vector<double>& v1,
                                      const vector<double>& p1,
                                      const vector<double>& v2,
                                      const vector<double>& p2,
                                      double d){
        
        double xx = (v1[0] + d*p1[0]) * (v2[0] + d*p2[0]);
        double yy = (v1[1] + d*p1[1]) * (v2[1] + d*p2[1]);
        double zz = (v1[2] + d*p1[2]) * (v2[2] + d*p2[2]);
        return xx + yy + zz;
        
    }
    
    
    
    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3). Returns a 3d vector.
    inline vector<double> vectorProduct(const vector<double>& v1,
                                        const vector<double>& v2,
                                        const vector<double>& v3,
                                        const vector<double>& v4) {
        vector<double> v;
        
        double vx = (v2[1]-v1[1])*(v4[2]-v3[2]) - (v2[2]-v1[2])*(v4[1]-v3[1]);
        double vy = (v2[2]-v1[2])*(v4[0]-v3[0]) - (v2[0]-v1[0])*(v4[2]-v3[2]);
        double vz = (v2[0]-v1[0])*(v4[1]-v3[1]) - (v2[1]-v1[1])*(v4[0]-v3[0]);
        
        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);
        
        return v;
    };
    
    
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z]. Returns a 3d vector.
    inline vector<double> crossProduct(const vector<double>& v1,
                                       const vector<double>& v2) {
        
        vector<double> v;
        
        double vx = v1[1]*v2[2] - v1[2]*v2[1];
        double vy = v1[2]*v2[0] - v1[0]*v2[2];
        double vz = v1[0]*v2[1] - v1[1]*v2[0];
        
        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);
        
        return v;
    };
    
    
    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3), but with v -> v+d*p. Returns a 3d vector.
    inline vector<double> vectorProductStretched (const vector<double>& v1,
                                                  const vector<double>& p1,
                                                  const vector<double>& v2,
                                                  const vector<double>& p2,
                                                  const vector<double>& v3,
                                                  const vector<double>& p3,
                                                  const vector<double>& v4,
                                                  const vector<double>& p4, double d){
        vector<double> v;
        
        double vx =
        ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]))
         - ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]));
        
        double vy =
        ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]))
        - ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]));
        
        double vz =
        ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]))
        - ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]));
        
        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);
        
        return v;
        
        
    };
    
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z]. Returns a 3d vector.
    inline vector<double> crossProductStretched(const vector<double>& v1,
                                                const vector<double>& p1,
                                                const vector<double>& v2,
                                                const vector<double>& p2,
                                                double d) {
        
        vector<double> v;
        
        double vx = (v1[1]+d*p1[1])*(v2[2]+d*p2[2]) - (v1[2]+d*p1[2])*(v2[1]+d*p2[1]);
        double vy = (v1[2]+d*p1[2])*(v2[0]+d*p2[0]) - (v1[0]+d*p1[0])*(v2[2]+d*p2[2]);
        double vz = (v1[0]+d*p1[0])*(v2[1]+d*p2[1]) - (v1[1]+d*p1[1])*(v2[0]+d*p2[0]);
        
        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);
        
        return v;
        
        
    };
    
    
    /// Projection of a new point based on a given direction and starting point
    inline vector<double> nextPointProjection(const vector<double>& coordinate,
                                              double d, const vector<double>& tau) {
        vector<double> v;
        v.push_back(coordinate[0] + d * tau[0]);
        v.push_back(coordinate[1] + d * tau[1]);
        v.push_back(coordinate[2] + d * tau[2]);
        return v;
    }
    
    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v| = alpha.
    inline vector<double> midPointCoordinate(const vector<double>& v1,
                                             const vector<double>& v2, double alpha) {
        vector<double> v;
        v.push_back(v1[0]*(1.0 - alpha) + alpha*v2[0]);
        v.push_back(v1[1]*(1.0 - alpha) + alpha*v2[1]);
        v.push_back(v1[2]*(1.0 - alpha) + alpha*v2[2]);
        return v;
    }
    
    
    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v| = alpha, but with x-d*p coordinates
    inline vector<double> midPointCoordinateStretched(const vector<double>& v1,
                                                      const vector<double>& p1,
                                                      const vector<double>& v2,
                                                      const vector<double>& p2,
                                                      double alpha , double d) {
        
        vector<double> v;
        v.push_back((v1[0] + d*p1[0])*(1.0 - alpha) + alpha*(v2[0] + d*p2[0]));
        v.push_back((v1[1] + d*p1[1])*(1.0 - alpha) + alpha*(v2[1] + d*p2[1]));
        v.push_back((v1[2] + d*p1[2])*(1.0 - alpha) + alpha*(v2[2] + d*p2[2]));
        return v;
    }
    
    /// Returns true if two vectors (p1->p2 and p3->p4) are parallel
    inline bool ifParallel(const vector<double>& p1, const vector<double>& p2,
                           const vector<double>& p3, const vector<double>& p4) {
        
        vector<double> v;
        
        v.push_back( (p2[2]-p1[2])*(p4[3]-p3[3]) -  (p2[3]-p1[3])*(p4[2]-p3[2]) );
        v.push_back( (p2[3]-p1[3])*(p4[1]-p3[1]) -  (p2[1]-p1[1])*(p4[3]-p3[3]) );
        v.push_back( (p2[1]-p1[1])*(p4[2]-p3[2]) -  (p2[2]-p1[2])*(p4[1]-p3[1]) );
        
        double norm = sqrt( v[1]*v[1] + v[2]*v[2] + v[3]*v[3] );
        return norm <= 1e-28 || norm == numeric_limits<double>::infinity();
    }
    
    /// Function to calculate a diatance between two segments
    double twoSegmentDistance(const vector<double>& v1, const vector<double>& v2,
                              const vector<double>& v3, const vector<double>& v4);
    
    /// Function to move bead out of plane by specified amount
    vector<double> movePointOutOfPlane(const vector<double>& p1,
                                       const vector<double>& p2,
                                       const vector<double>& p3,
                                       const vector<double>& p4, int i, double d);
    
    
    /// Function to create a initial branching point and direction, given an
    /// initial normal vector and point.
    /// @param l - the distance of the branch from the original point
    /// @param m - the size of the branch projection
    /// @param theta - the angle of branching
    /// @return a vector describing the initial branching direction and point
    tuple<vector<double>, vector<double>> branchProjection(const vector<double>& n,
                                                           const vector<double>& p,
                                                           double l, double m, double theta);
    
    
    }

#endif
