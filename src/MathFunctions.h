
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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

// A simple coordinate type that makes operations easier
// Can be replaced by linalg/tensor libraries in the future
template< size_t dim, typename Float = double > struct Vec {
    std::array< Float, dim > value;

    using storage_type = std::array< Float, dim >;
    using size_type = typename storage_type::size_type;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    constexpr iterator       begin()       noexcept { return value.begin(); }
    constexpr const_iterator begin() const noexcept { return value.begin(); }
    constexpr iterator       end()       noexcept { return value.end(); }
    constexpr const_iterator end() const noexcept { return value.end(); }

    constexpr       Float& operator[](size_type pos)       { return value[pos]; }
    constexpr const Float& operator[](size_type pos) const { return value[pos]; }
};

using Vec3 = Vec<3>;
    
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

template< size_t dim, typename Float >
inline Float magnitude2(const Vec<dim, Float>& v) {
    Float mag2 = 0.0;
    for(size_t i = 0; i < dim; ++i) mag2 += v[i] * v[i];
    return mag2;
}
template< size_t dim, typename Float >
inline Float magnitude(const Vec<dim, Float>& v) {
    return sqrt(magnitude2(v));
}

template< size_t dim, typename Float >
inline void normalize(Vec<dim, Float>& v) {
    Float norm = magnitude(v);
    for(size_t i = 0; i < dim; ++i) v[i] /= norm;
}
template< size_t dim, typename Float >
inline auto normalizedVector(const Vec<dim, Float>& v) {
    Vec<dim, Float> res = v;
    normalize(res);
    return res;
}
    
    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3)
    inline double twoPointDistance(const vector<double>& v1, const vector<double>& v2) {
        
        return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) +
                    (v2[1]-v1[1])*(v2[1]-v1[1]) +
                    (v2[2]-v1[2])*(v2[2]-v1[2]));
    }
template< size_t dim, typename Float >
inline Float distance2(const Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    Float res = 0.0;
    for(size_t idx = 0; idx < dim; ++idx) {
        res += (v2[idx] - v1[idx]) * (v2[idx] - v1[idx]);
    }
    return res;
}
template< size_t dim, typename Float >
inline Float distance(const Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    return sqrt(distance2(v1, v2));
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
template< size_t dim, typename Float >
inline Float dot(const Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    Float res = 0.0;
    for(size_t idx = 0; idx < dim; ++idx)
        res += v1[idx] * v2[idx];
    return res;
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

    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    inline double scalarProduct(const vector<double>& v1, const vector<double>& v2,
                                const vector<double>& v3, const vector<double>& v4) {
        
        return ((v2[0] - v1[0])*(v4[0] - v3[0]) +
                (v2[1] - v1[1])*(v4[1] - v3[1]) +
                (v2[2] - v1[2])*(v4[2] - v3[2]));
    }
    template< size_t Dim >
    inline double scalarProduct(
        const array<double, Dim>& v1, const array<double, Dim>& v2,
        const array<double, Dim>& v3, const array<double, Dim>& v4
    ) {
        double ret = 0.0;
        for(size_t i = 0; i < Dim; ++i) {
            ret += (v2[i] - v1[i]) * (v4[i] - v3[i]);
        }
        return ret;
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
    }
    inline array<double, 3> vectorProduct(const array<double, 3>& v1,
                                          const array<double, 3>& v2,
                                          const array<double, 3>& v3,
                                          const array<double, 3>& v4) {
        return array<double, 3> {{
            (v2[1]-v1[1])*(v4[2]-v3[2]) - (v2[2]-v1[2])*(v4[1]-v3[1]),
            (v2[2]-v1[2])*(v4[0]-v3[0]) - (v2[0]-v1[0])*(v4[2]-v3[2]),
            (v2[0]-v1[0])*(v4[1]-v3[1]) - (v2[1]-v1[1])*(v4[0]-v3[0])
        }};
    }
    
    
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
inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
    return Vec3 {{{
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    }}};
}
    
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
template< size_t dim, typename Float >
inline auto vector2Vec(const vector<Float>& v) {
    // Assert v.size() == Size
    Vec<dim, Float> res;
    for(size_t idx = 0; idx < dim; ++idx){
        res[idx] = v[idx];
    }
    return res;
}
template< size_t dim, typename Float >
inline auto vec2Vector(const Vec<dim, Float>& a) {
    return vector<Float>(a.value.begin(), a.value.end());
}

template< size_t dim, typename Float >
inline auto operator-(const Vec<dim, Float>& v){
    Vec<dim, Float> res;
    for(size_t idx = 0; idx < dim; ++idx){
        res[idx] = -v[idx];
    }
    return res;
}
template< size_t dim, typename Float >
inline auto operator+(const Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    Vec<dim, Float> res;
    for(size_t idx1 = 0; idx1 < dim; ++idx1) {
        res[idx1] = v1[idx1] + v2[idx1];
    }
    return res;
}
template< size_t dim, typename Float >
inline auto& operator+=(Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    for(size_t idx = 0; idx < dim; ++idx) {
        v1[idx] += v2[idx];
    }
    return v1;
}
template< size_t dim, typename Float >
inline auto operator-(const Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    Vec<dim, Float> res;
    for(size_t idx1 = 0; idx1 < dim; ++idx1) {
        res[idx1] = v1[idx1] - v2[idx1];
    }
    return res;
}
template< size_t dim, typename Float >
inline auto& operator-=(Vec<dim, Float>& v1, const Vec<dim, Float>& v2) {
    for(size_t idx = 0; idx < dim; ++idx) {
        v1[idx] -= v2[idx];
    }
    return v1;
}
template< size_t dim, typename Float >
inline auto operator*(const Vec<dim, Float>& v, Float k) {
    Vec<dim, Float> res;
    for(size_t idx = 0; idx < dim; ++idx){
        res[idx] = v[idx] * k;
    }
    return res;
}
template< size_t dim, typename Float >
inline auto operator*(Float k, const Vec<dim, Float>& v) {
    return v * k;
}
template< size_t dim, typename Float >
inline auto& operator*=(Vec<dim, Float>& v, double k) {
    for(size_t i = 0l i < dim; ++i) v[i] *= k;
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
    constexpr array<array<double, 3>, 3> Eye3 = {{ {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}} }};

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
