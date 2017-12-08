
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_MathFunctions_h
#define MEDYAN_MathFunctions_h

#include <array>
#include <cmath>
#include <vector>

#include "common.h"

/// @namespace mathfunc is used for the mathematics module for the entire codebase
/// mathfunc includes functions to calculate distances, products, and midpoints

namespace mathfunc {
    
    /// Normalize a vector
    inline void normalize(vector<double>& v) {
        
        double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        
        v[0] /= norm; v[1] /= norm; v[2] /= norm;
    }
    
    /// Return normalized vector
    inline vector<double> normalizedVector(const vector<double>& v) {
        
        double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        
        vector<double> v1;
        
        v1.push_back(v[0]/norm); v1.push_back(v[1]/norm); v1.push_back(v[2]/norm);
        
        return v1;
    }
    
    /// Get the magnitude of a vector
    inline double magnitude(const vector<double>& v) {
        
        return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    template<size_t Dim>
    inline double magnitude(const array<double, Dim>& v) {
        double mag2 = 0;
        for(auto it: v) mag2 += it * it;
        return sqrt(mag2);
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
    
    /// Scalar product of two vectors v1(x,y,z) and v2(x,y,z)
    inline double dotProduct(const vector<double>& v1, const vector<double>& v2) {
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }
    
    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    inline double scalarProduct(const vector<double>& v1, const vector<double>& v2,
                                const vector<double>& v3, const vector<double>& v4) {
        
        return ((v2[0] - v1[0])*(v4[0] - v3[0]) +
                (v2[1] - v1[1])*(v4[1] - v3[1]) +
                (v2[2] - v1[2])*(v4[2] - v3[2]));
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
                                         const vector<double>& p4,
                                         double d){
        
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
    
    
    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3), but with v -> v+d*p. Returns a 3d vector.
    inline vector<double> vectorProductStretched (const vector<double>& v1,
                                                  const vector<double>& p1,
                                                  const vector<double>& v2,
                                                  const vector<double>& p2,
                                                  const vector<double>& v3,
                                                  const vector<double>& p3,
                                                  const vector<double>& v4,
                                                  const vector<double>& p4,
                                                  double d){
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

    /// Vector and array converter. Need to ensure the vector has size of _Size
    // No need for move semantics because normally we use this for copying integers or doubles
    template<typename Ty, size_t Size>
    inline array<Ty, Size> vector2Array(const vector<Ty>& v) {
        // Assert v.size() == Size
        array<Ty, Size> res;
        for(size_t idx = 0; idx < Size; ++idx){
            res[idx] = v[idx];
        }
        return res;
    }
    template<typename Ty, size_t Size>
    inline vector<Ty> array2Vector(const array<Ty, Size>& a) {
        return vector<Ty>(a.begin(), a.end());
    }

    /// Get the negative of the vector
    template<size_t Dim>
    inline array<double, Dim> vectorNegative(const array<double, Dim>& v){
        array<double, Dim> res;
        for(size_t idx = 0; idx < Dim; ++idx){
            res[idx] = -v[idx];
        }
        return res;
    }
    /// Vector sum.
    template<size_t Dim>
    inline array<double, Dim> vectorSum(const array<double, Dim>& v1,
                                        const array<double, Dim>& v2) {
        array<double, Dim> res;
        for(size_t idx1 = 0; idx1 < Dim; ++idx1) {
            res[idx1] = v1[idx1] + v2[idx1];
        }
        return res;
    }
    /// Vector difference.
    template<size_t Dim>
    inline array<double, Dim> vectorDifference(const array<double, Dim>& v1,
                                               const array<double, Dim>& v2) {
        array<double, Dim> res;
        for(size_t idx1 = 0; idx1 < Dim; ++idx1) {
            res[idx1] = v1[idx1] - v2[idx1];
        }
        return res;
    }
    /// Vector multiply.
    template<size_t Dim>
    inline array<double, Dim> vectorMultiply(const array<double, Dim>& v,
                                             const double k) {
        array<double, Dim> res;
        for(size_t idx = 0; idx < Dim; ++idx){
            res[idx] = v[idx] * k;
        }
        return res;
    }
    /// Vector expand
    template<size_t Dim>
    inline array<double, Dim>& vectorExpand(array<double, Dim>& v,
                                            const double k) {
        for(auto& it: v) it *= k;
        return v;
    }

    /// Get the negative of the matrix
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<array<double, Dim2>, Dim1> matrixNegative(const array<array<double, Dim2>, Dim1>& m){
        array<array<double, Dim2>, Dim1> res;

        for(size_t idx1 = 0; idx1 < Dim1; ++idx1){
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2){
                res[idx1][idx2] = -m[idx1][idx2];
            }
        }
        return res;
    }
    /// Matrix sum. MUST have same dimension
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<array<double, Dim2>, Dim1> matrixSum(const array<array<double, Dim2>, Dim1>& m1,
                                                      const array<array<double, Dim2>, Dim1>& m2) {
        array<array<double, Dim2>, Dim1> res;
        for(size_t idx1 = 0; idx1 < Dim1; ++idx1){
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2){
                res[idx1][idx2] = m1[idx1][idx2] + m2[idx1][idx2];
            }
        }
        return res;
    }
    /// Add matrix2 to matrix1. Returns a reference to the modified vector
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<array<double, Dim2>, Dim1>& matrixIncrease(array<array<double, Dim2>, Dim1>& m1,
                                                      const array<array<double, Dim2>, Dim1>& m2) {
        for(size_t idx1 = 0; idx1 < Dim1; ++idx1){
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2){
                m1[idx1][idx2] += m2[idx1][idx2];
            }
        }
        return m1;
    }
    /// Matrix difference. MUST have same dimension
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<array<double, Dim2>, Dim1> matrixDifference(const array<array<double, Dim2>, Dim1>& m1,
                                                             const array<array<double, Dim2>, Dim1>& m2) {
        array<array<double, Dim2>, Dim1> res;
        for(size_t idx1 = 0; idx1 < Dim1; ++idx1){
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2){
                res[idx1][idx2] = m1[idx1][idx2] - m2[idx1][idx2];
            }
        }
        return res;
    }
    /// Matrix multiply.
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<array<double, Dim2>, Dim1> matrixMultiply(const array<array<double, Dim2>, Dim1>& m,
                                                           const double k) {
        array<array<double, Dim2>, Dim1> res;
        for(size_t idx1 = 0; idx1 < Dim1; ++idx1){
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2){
                res[idx1][idx2] = m[idx1][idx2] * k;
            }
        }
        return res;
    }

    /// Identity matrix 3x3. Internal linkage applied implicitly by using "const".
    const array<array<double, 3>, 3> Eye3 = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    /// Tensor product of two vectors, notated as a matrix with size dim(v1) x dim(v2)
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<array<double, Dim2>, Dim1> tensorProduct(const array<double, Dim1>& v1,
                                                          const array<double, Dim2>& v2) {
        array<array<double, Dim2>, Dim1> res;
        for(size_t idx1 = 0; idx1 < Dim1; ++idx1) {
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2) {
                res[idx1][idx2] = v1[idx1] * v2[idx2];
            }
        }
        return res;
    }

    /// Matrix (d1*d2) times vector (d2). Returns a vector with dimension d1
    template<size_t Dim1, size_t Dim2=Dim1>
    inline array<double, Dim1> matrixProduct(const array<array<double, Dim2>, Dim1>& m,
                                             const array<double, Dim2>& v) {
        array<double, Dim1> res = {};
        
        for(size_t idx1 = 0; idx1 < Dim1; ++idx1){
            for(size_t idx2 = 0; idx2 < Dim2; ++idx2){
                res[idx1] += m[idx1][idx2] * v[idx2];
            }
        }

        return res;
    }

    /// Returns the area of a triangle with vertices at point v1, v2 and v3.
    inline double areaTriangle(const vector<double>& v1,
                               const vector<double>& v2,
                               const vector<double>& v3) {
        
        return 0.5 * magnitude(vectorProduct(v1, v2, v1, v3));
    }

    /// The area of a triangle with vertices at point v1 + d*p1, v2 + d*p2 and v3 + d*p3
    inline double areaTriangleStretched(const vector<double>& v1,
                                        const vector<double>& p1,
                                        const vector<double>& v2,
                                        const vector<double>& p2,
                                        const vector<double>& v3,
                                        const vector<double>& p3,
                                        double d) {
        return 0.5 * magnitude(vectorProductStretched(v1, p1, v2, p2, v1, p1, v3, p3, d));

    }

    /// The gradient on the area of a (triangle with vertices at point v1, v2 and v3)
    /// with regard to coordinate v1, v2 and v3
    inline array<array<double, 3>, 3> areaTriangleGradient(const vector<double>& v1,
                                                           const vector<double>& v2,
                                                           const vector<double>& v3) {
        double area = areaTriangle(v1, v2, v3);

        array<array<double, 3>, 3> result;

        double dot23 = scalarProduct(v1, v2, v1, v3);
        double r12 = twoPointDistance(v1, v2), r23 = twoPointDistance(v1, v3);

        result[0][0] = (-r12*r12*(v3[0] - v1[0]) - r23*r23*(v2[0] - v1[0]) + dot23*(v2[0] + v3[0] - 2*v1[0])) / area / 4;
        result[0][1] = (-r12*r12*(v3[1] - v1[1]) - r23*r23*(v2[1] - v1[1]) + dot23*(v2[1] + v3[1] - 2*v1[1])) / area / 4;
        result[0][2] = (-r12*r12*(v3[2] - v1[2]) - r23*r23*(v2[2] - v1[2]) + dot23*(v2[2] + v3[2] - 2*v1[2])) / area / 4;
        
        result[1][0] = (r23*r23*(v2[0] - v1[0]) - dot23*(v3[0] - v1[0])) / area / 4;
        result[1][1] = (r23*r23*(v2[1] - v1[1]) - dot23*(v3[1] - v1[1])) / area / 4;
        result[1][2] = (r23*r23*(v2[2] - v1[2]) - dot23*(v3[2] - v1[2])) / area / 4;

        result[2][0] = (r12*r12*(v3[0] - v1[0]) - dot23*(v2[0] - v1[0])) / area / 4;
        result[2][1] = (r12*r12*(v3[1] - v1[1]) - dot23*(v2[1] - v1[1])) / area / 4;
        result[2][2] = (r12*r12*(v3[2] - v1[2]) - dot23*(v2[2] - v1[2])) / area / 4;        

        return result;
        // TODO: Test the result!!!!!!!!
    }
    
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
    /// |v-v1|/|v2-v1| = alpha.
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
    inline bool areParallel(const vector<double>& p1, const vector<double>& p2,
                            const vector<double>& p3, const vector<double>& p4) {
        
        auto v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        auto v2 = {p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]};

        return areEqual(magnitude(crossProduct(v1,v2)), 0.0);
    }
    
    /// Returns true if two vectors (p1->p2 and p3->p4) are in the same plane
    inline bool areInPlane(const vector<double>& p1, const vector<double>& p2,
                           const vector<double>& p3, const vector<double>& p4) {
        
        auto v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        auto v2 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};
        auto v3 = {p4[0] - p1[0], p4[1] - p1[1], p4[2] - p1[2]};
        
        auto cp = crossProduct(v1, v2);
        
        return areEqual(dotProduct(v3, cp), 0.0);
    }
    
    /// Function to move bead out of plane by specified amount
    vector<double> movePointOutOfPlane(const vector<double>& p1,
                                       const vector<double>& p2,
                                       const vector<double>& p3,
                                       const vector<double>& p4,
                                       int i, double d);
    
    
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
