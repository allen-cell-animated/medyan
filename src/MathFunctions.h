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

#include <cmath>
#include <vector>

#include "common.h"
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#endif

/// @namespace mathfunc is used for the mathematics module for the entire codebase
/// mathfunc includes functions to calculate distances, products, and midpoints

namespace mathfunc {
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

#else
    static __inline__ __device__ double atomicAdd(double *address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    if (val==0.0)
      return __longlong_as_double(old);
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
  }
    static __inline__ __device__ double atomicExch(double *address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    if (val==0.0)
      return __longlong_as_double(old);
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(__longlong_as_double(0.0)));
    } while (assumed != old);
    return __longlong_as_double(old);
  }
#endif
#ifdef CUDAACCL
     __host__ inline int nextPowerOf2( int n)
    {
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;
        return n;
    }
    __host__ inline vector<int> getaddred2bnt(int nint){
        vector<int> bnt;
        int THREADSPERBLOCK = 0;
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        THREADSPERBLOCK = prop.maxThreadsPerBlock;
        auto maxthreads = 8 * THREADSPERBLOCK;
        int M = max(nextPowerOf2(nint),64*4);
        int blocks, threads;
        if(M > THREADSPERBLOCK){
            if(M > maxthreads) {
                blocks = 8;
                threads = THREADSPERBLOCK;
            }
            else if(M > THREADSPERBLOCK){
                blocks = M /(4 * THREADSPERBLOCK) +1;
                threads = THREADSPERBLOCK;
            }
        }
        else
        { blocks = 1; threads = M/4;}
        //0 M
        //1 THREADSPERBLOCK
        //2 blocks
        //3 threads
        bnt.clear();
        bnt.push_back(M);
        bnt.push_back(THREADSPERBLOCK);
        bnt.push_back(blocks);
        bnt.push_back(threads);
        return bnt;
    }
    __global__ void resetintvariableCUDA(int *variable);
    __global__ void resetdoublevariableCUDA(double *variable);
     __global__ void addvector(double *U, int *params, double *U_sum, double *U_tot);
//    __global__ void addvector(double *U, int *params, double *U_sum, double *U_tot, int* culpritID, char* culpritFF,
//                              char* culpritinteraction, char *FF, char *interaction);
//    __global__ void addvectorred(double *U, int *params, double *U_sum, double *U_tot);
    __global__ void addvectorred2(double *U, int *params, double *U_sum, double *U_tot);
    #endif
    /// Normalize a vector
    inline void normalize(vector<double> &v) {

        double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }

    /// Return normalized vector not in place
    inline vector<double> normalizeVector(const vector<double> &v) {

        double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

        vector<double> v1;

        v1.push_back(v[0] / norm);
        v1.push_back(v[1] / norm);
        v1.push_back(v[2] / norm);

        return v1;
    }

    /// Return normalized vector
    /// ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void normalizeVector(double *v) {

        double norm = sqrt((*(v)) * (*(v)) + (*(v + 1)) * (*(v + 1)) + (*(v + 2)) * (*(v + 2)));
        *(v) = *(v) / norm;
        *(v + 1) = *(v + 1) / norm;
        *(v + 2) = *(v + 2) / norm;
    }

    /// Get the magnitude of a vector
    inline double magnitude(const vector<double> &v) {

        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    ///ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double magnitude(double const *v) {

        return sqrt((*(v)) * (*(v)) + (*(v + 1)) * (*(v + 1)) + (*(v + 2)) * (*(v + 2)));
    }

    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double getdistancefromplane(double *coord, double * plane,int id){
        return (plane[0] * coord[id] + plane[1] * coord[id + 1] + plane[2] * coord[id + 2] + plane[3]) /
               sqrt(pow(plane[0], 2) + pow(plane[1], 2) + pow(plane[2], 2));
    }

    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double getstretcheddistancefromplane(double *coord, double *force, double * plane, double z, int id){
        double movedPoint[3] = {coord[id] + z*force[id],
                                     coord[id + 1] + z*force[id + 1],
                                     coord[id + 2] + z*force[id + 2]};
        return (plane[0] * movedPoint[0] + plane[1] * movedPoint[1] + plane[2] * movedPoint[2] + plane[3]) /
               sqrt(pow(plane[0], 2) + pow(plane[1], 2) + pow(plane[2], 2));
    }
    //@{
    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3)
    inline double twoPointDistance(const vector<double> &v1, const vector<double> &v2) {

        return sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                    (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    inline double twoPointDistance(double const *v1, const vector<double> &v2) {

        return sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                    (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    inline double twoPointDistance(const vector<double> &v1, double const *v2) {

        return sqrt((*(v2) - v1[0]) * (*(v2) - v1[0]) +
                    (*(v2 + 1) - v1[1]) * (*(v2 + 1) - v1[1]) +
                    (*(v2 + 2) - v1[2]) * (*(v2 + 2) - v1[2]));
    }
    //@}

    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3)
    /// ARRAY VERSION
    #ifdef CUDAACCL
     __host__ __device__
     #endif
    inline double twoPointDistance(double const *v1, double const *v2) {

        return sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double twoPointDistance(double const *v1, double const *v2, int const id) {

        return sqrt((v2[id] - v1[id]) * (v2[id] - v1[id]) + (v2[id + 1] - v1[id + 1]) * (v2[id + 1] - v1[id + 1])
                    + (v2[id + 2] - v1[id + 2]) * (v2[id + 2] - v1[id + 2]));
    }
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double twoPointDistancemixedID(double const *v1, double const *v2, int const id1, int const id2) {

        return sqrt((v2[id2] - v1[id1]) * (v2[id2] - v1[id1]) + (v2[id2 + 1] - v1[id1 + 1]) * (v2[id2 + 1] - v1[id1 +
                1]) + (v2[id2 + 2] - v1[id1 + 2]) * (v2[id2 + 2] - v1[id1 + 2]));
    }
    /// Compute distance between two points with coordinates
    /// (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z)
    inline double twoPointDistanceStretched(const vector<double> &v1,
                                            const vector<double> &p1,
                                            const vector<double> &v2,
                                            const vector<double> &p2, double d) {

//        std::cout<<p1[0]<<" "<<p1[1]<<" "<<p1[2]<<" "<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<<endl;
//        std::cout<<((v2[0] + d*p2[0])-(v1[0] + d*p1[0]))<<" "<<((v2[0] + d*p2[0])-(v1[0] + d*p1[0]))<<" "<<((v2[1] + d*p2[1])-(v1[1] + d*p1[1]))<<
//        " "<<((v2[1] + d*p2[1])-(v1[1] + d*p1[1]))<<" "<<((v2[2] + d*p2[2])-(v1[2] + d*p1[2]))<<" "<<((v2[2] + d*p2[2])-(v1[2] + d*p1[2]))<<endl;
        return sqrt(((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) +
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) +
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])));
    }

    /// Compute distance between two points with coordinates
    /// (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z)
    /// ARRAY VERSION
    inline double twoPointDistanceStretched(double const *v1,
                                            double const *p1,
                                            double const *v2,
                                            double const *p2, double d) {
//                std::cout<<p1[0]<<" "<<p1[1]<<" "<<p1[2]<<" "<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<<endl;
//        std::cout<<((v2[0] + d*p2[0])-(v1[0] + d*p1[0]))<<" "<<((v2[0] + d*p2[0])-(v1[0] + d*p1[0]))<<" "<<((v2[1] + d*p2[1])-(v1[1] + d*p1[1]))<<
//        " "<<((v2[1] + d*p2[1])-(v1[1] + d*p1[1]))<<" "<<((v2[2] + d*p2[2])-(v1[2] + d*p1[2]))<<" "<<((v2[2] + d*p2[2])-(v1[2] + d*p1[2]))<<endl;
        return sqrt(((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) +
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) +
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])));
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double twoPointDistanceStretched(double const *v1,
                                            double const *p1,
                                            double const *v2,
                                            double const *p2, double d, int const id) {
        return sqrt(((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                    ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) +
                    ((v2[id + 1] + d * p2[id + 1]) - (v1[id + 1] + d * p1[id + 1])) *
                    ((v2[id + 1] + d * p2[id + 1]) - (v1[id + 1] + d * p1[id + 1])) +
                    ((v2[id + 2] + d * p2[id + 2]) - (v1[id + 2] + d * p1[id + 2])) *
                    ((v2[id + 2] + d * p2[id + 2]) - (v1[id + 2] + d * p1[id + 2])));
    }
    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double twoPointDistanceStretchedmixedID(double const *v1,
                                            double const *p1,
                                            double const *v2,
                                            double const *p2, double d, int const id1, const int id2) {
        return sqrt(((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                    ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) +
                    ((v2[id2 + 1] + d * p2[id2 + 1]) - (v1[id1 + 1] + d * p1[id1 + 1])) *
                    ((v2[id2 + 1] + d * p2[id2 + 1]) - (v1[id1 + 1] + d * p1[id1 + 1])) +
                    ((v2[id2 + 2] + d * p2[id2 + 2]) - (v1[id1 + 2] + d * p1[id1 + 2])) *
                    ((v2[id2 + 2] + d * p2[id2 + 2]) - (v1[id1 + 2] + d * p1[id1 + 2])));
    }
    //@{
    /// Calculates a normal to a line starting at (x1,y1,z1) and ending at (x2,y2,z2)
    inline vector<double> twoPointDirection(const vector<double> &v1,
                                            const vector<double> &v2) {
        vector<double> tau(3, 0);
        double invD = 1 / twoPointDistance(v1, v2);
        tau[0] = invD * (v2[0] - v1[0]);
        tau[1] = invD * (v2[1] - v1[1]);
        tau[2] = invD * (v2[2] - v1[2]);
        return tau;
    }

    inline vector<double> twoPointDirection(double const *v1,
                                            const vector<double> &v2) {
        vector<double> tau(3, 0);
        double invD = 1 / twoPointDistance(v1, v2);
        tau[0] = invD * (v2[0] - v1[0]);
        tau[1] = invD * (v2[1] - v1[1]);
        tau[2] = invD * (v2[2] - v1[2]);
        return tau;
    }

    inline void twoPointDirection(double *tau,
                                  double const *v1,
                                  double const *v2) {

        double invD = 1 / twoPointDistance(v1, v2);
        tau[0] = invD * (*(v2) - *(v1));
        tau[1] = invD * (*(v2 + 1) - *(v1 + 1));
        tau[2] = invD * (*(v2 + 2) - *(v1 + 2));
    }
    //@}

    /// Scalar product of two vectors v1(x,y,z) and v2(x,y,z)
    inline double dotProduct(const vector<double> &v1, const vector<double> &v2) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    /// Scalar product of two vectors v1(x,y,z) and v2(x,y,z)
    /// ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double dotProduct(double const *v1, double const *v2) {
        return (*(v1)) * (*(v2)) + (*(v1 + 1)) * (*(v2 + 1)) + (*(v1 + 2)) * (*(v2 + 2));

//        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }

    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    inline double scalarProduct(const vector<double> &v1, const vector<double> &v2,
                                const vector<double> &v3, const vector<double> &v4) {

        return ((v2[0] - v1[0]) * (v4[0] - v3[0]) +
                (v2[1] - v1[1]) * (v4[1] - v3[1]) +
                (v2[2] - v1[2]) * (v4[2] - v3[2]));
    }

    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    /// ARRAY VERSION
    inline double scalarProduct(double const *v1, double const *v2,
                                double const *v3, double const *v4) {
//        return((*(v2)-*(v1))*(*(v4)-*(v3))
//               +(*(v2+1)-*(v1+1))*(*(v4+1)-*(v3+1))
//               +(*(v2+2)-*(v1+2))*(*(v4+2)-*(v3+2)));
        return ((v2[0] - v1[0]) * (v4[0] - v3[0]) +
                (v2[1] - v1[1]) * (v4[1] - v3[1]) +
                (v2[2] - v1[2]) * (v4[2] - v3[2]));
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double scalarProduct(double const *v1, double const *v2,
                                double const *v3, double const *v4, int const id) {
//        return((*(v2)-*(v1))*(*(v4)-*(v3))
//               +(*(v2+1)-*(v1+1))*(*(v4+1)-*(v3+1))
//               +(*(v2+2)-*(v1+2))*(*(v4+2)-*(v3+2)));

        return ((v2[id] - v1[id]) * (v4[id] - v3[id]) +
                (v2[id + 1] - v1[id + 1]) * (v4[id + 1] - v3[id + 1]) +
                (v2[id + 2] - v1[id + 2]) * (v4[id + 2] - v3[id + 2]));
    }
    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline double scalarProductmixedID(double const *v1, double const *v2,
                                double const *v3, double const *v4, int const id1, int const id2, int const id3, int
                                       const id4) {
        return ((v2[id2] - v1[id1]) * (v4[id4] - v3[id3]) +
                (v2[id2 + 1] - v1[id1 + 1]) * (v4[id4 + 1] - v3[id3 + 1]) +
                (v2[id2 + 2] - v1[id1 + 2]) * (v4[id4 + 2] - v3[id3 + 2]));
    }
    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3) but with x+d*p coordinates
    inline double scalarProductStretched(const vector<double> &v1,
                                         const vector<double> &p1,
                                         const vector<double> &v2,
                                         const vector<double> &p2,
                                         const vector<double> &v3,
                                         const vector<double> &p3,
                                         const vector<double> &v4,
                                         const vector<double> &p4,
                                         double d) {

        double xx = ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v4[0] + d * p4[0]) - (v3[0] + d * p3[0]));
        double yy = ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v4[1] + d * p4[1]) - (v3[1] + d * p3[1]));
        double zz = ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v4[2] + d * p4[2]) - (v3[2] + d * p3[2]));
        return xx + yy + zz;

    }

    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3) but with x+d*p coordinates
    ///
    inline double scalarProductStretched(double const *v1,
                                         double const *p1,
                                         double const *v2,
                                         double const *p2,
                                         double const *v3,
                                         double const *p3,
                                         double const *v4,
                                         double const *p4,
                                         double d) {
        double xx = ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v4[0] + d * p4[0]) - (v3[0] + d * p3[0]));
        double yy = ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v4[1] + d * p4[1]) - (v3[1] + d * p3[1]));
        double zz = ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v4[2] + d * p4[2]) - (v3[2] + d * p3[2]));
        return xx + yy + zz;
//        double xx = ((*(v2) + d*(*(p2))) - (*(v1) + d*(*(p1)))) *
//        ((*(v4)) + d*(*(p4))) - ((*(v3)) + d*(*(p3)));
//        double yy = ((*(v2+1) + d*(*(p2+1))) - (*(v1+1) + d*(*(p1+1))))*
//        ((*(v4+1) + d*(*(p4+1))) - (*(v3+1) + d*(*(p3+1))));
//        double zz = ((*(v2+2) + d*(*(p2+2))) - (*(v1+2) + d*(*(p1+2))))*
//        ((*(v4+2) + d*(*(p4+2))) - (*(v3+2) + d*(*(p3+2))));
//        return xx + yy + zz;

    }
    //CUDA version
    #ifdef CUDAACCL
     __host__ __device__
     #endif
    inline double scalarProductStretched(double const *v1,
                                         double const *p1,
                                         double const *v2,
                                         double const *p2,
                                         double const *v3,
                                         double const *p3,
                                         double const *v4,
                                         double const *p4,
                                         double d,int const id ) {
        double xx = ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                    ((v4[id] + d * p4[id]) - (v3[id] + d * p3[id]));
        double yy = ((v2[id + 1] + d * p2[id + 1]) - (v1[id + 1] + d * p1[id + 1])) *
                    ((v4[id + 1] + d * p4[id + 1]) - (v3[id + 1] + d * p3[id + 1]));
        double zz = ((v2[id + 2] + d * p2[id + 2]) - (v1[id + 2] + d * p1[id + 2])) *
                    ((v4[id + 2] + d * p4[id + 2]) - (v3[id + 2] + d * p3[id + 2]));
        return xx + yy + zz;
    }

    #ifdef CUDAACCL
    __host__ __device__
     #endif
    inline double scalarProductStretchedmixedIDv2(double const *v1,

                                         double const *v2,

                                         double const *v3,

                                         double const *v4,

                                         double d, int const id, double const *p, int const id1, int const id2, int
                                                  const id3, int const id4) {
        double xx = ((v2[id] + d * p[id2]) - (v1[id] + d * p[id1])) *
                    ((v4[id] + d * p[id4]) - (v3[id] + d * p[id3]));
        double yy = ((v2[id + 1] + d * p[id2 + 1]) - (v1[id + 1] + d * p[id1 + 1])) *
                    ((v4[id + 1] + d * p[id4 + 1]) - (v3[id + 1] + d * p[id3 + 1]));
        double zz = ((v2[id + 2] + d * p[id2 + 2]) - (v1[id + 2] + d * p[id1 + 2])) *
                    ((v4[id + 2] + d * p[id4 + 2]) - (v3[id + 2] + d * p[id3 + 2]));
        return xx + yy + zz;
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
     #endif
    inline double scalarProductStretchedmixedID(double const *v1,
                                         double const *p1,
                                         double const *v2,
                                         double const *p2,
                                         double const *v3,
                                         double const *p3,
                                         double const *v4,
                                         double const *p4,
                                         double d,
                                         int const id1,
                                         int const id2,
                                         int const id3,
                                         int const id4) {
        double xx = ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                    ((v4[id4] + d * p4[id4]) - (v3[id3] + d * p3[id3]));
        double yy = ((v2[id2 + 1] + d * p2[id2 + 1]) - (v1[id1 + 1] + d * p1[id1 + 1])) *
                    ((v4[id4 + 1] + d * p4[id4 + 1]) - (v3[id3 + 1] + d * p3[id3 + 1]));
        double zz = ((v2[id2 + 2] + d * p2[id2 + 2]) - (v1[id1 + 2] + d * p1[id1 + 2])) *
                    ((v4[id4 + 2] + d * p4[id4 + 2]) - (v3[id3 + 2] + d * p3[id3 + 2]));
        return xx + yy + zz;
    }
    /// Scalar product of two vectors with coordinates: v1[z,y,z] + d*p1[x,y,z] and
    /// v2[x,y,z] + d*p2[x,y,z]
    inline double dotProductStretched(const vector<double> &v1,
                                      const vector<double> &p1,
                                      const vector<double> &v2,
                                      const vector<double> &p2,
                                      double d) {

        double xx = (v1[0] + d * p1[0]) * (v2[0] + d * p2[0]);
        double yy = (v1[1] + d * p1[1]) * (v2[1] + d * p2[1]);
        double zz = (v1[2] + d * p1[2]) * (v2[2] + d * p2[2]);
        return xx + yy + zz;

    }


    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3). Returns a 3d vector.
    inline vector<double> vectorProduct(const vector<double> &v1,
                                        const vector<double> &v2,
                                        const vector<double> &v3,
                                        const vector<double> &v4) {
        vector<double> v;

        double vx = (v2[1] - v1[1]) * (v4[2] - v3[2]) - (v2[2] - v1[2]) * (v4[1] - v3[1]);
        double vy = (v2[2] - v1[2]) * (v4[0] - v3[0]) - (v2[0] - v1[0]) * (v4[2] - v3[2]);
        double vz = (v2[0] - v1[0]) * (v4[1] - v3[1]) - (v2[1] - v1[1]) * (v4[0] - v3[0]);

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;
    };

    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3). Returns a 3d vector.
    /// ARRAY VERSION
    inline void vectorProduct(double *v,
                              double const *v1,
                              double const *v2,
                              double const *v3,
                              double const *v4) {
        double vx =
                (*(v2 + 1) - *(v1 + 1)) * (*(v4 + 2) - *(v3 + 2)) - (*(v2 + 2) - *(v1 + 2)) * (*(v4 + 1) - *(v3 + 1));
        double vy = (*(v2 + 2) - *(v1 + 2)) * (*(v4) - *(v3)) - (*(v2) - *(v1)) * (*(v4 + 2) - *(v3 + 2));
        double vz = (*(v2) - *(v1)) * (*(v4 + 1) - *(v3 + 1)) - (*(v2 + 1) - *(v1 + 1)) * (*(v4) - *(v3));

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };

    /// CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProduct(double *v,
                              double const *v1,
                              double const *v2,
                              double const *v3,
                              double const *v4, int const id) {
        double vx =
                (v2[id+1] - v1[id+1]) * (v4[id+2] - v3[id+2]) - (v2[id+2] - v1[id+2]) * (v4[id+1] - v3[id+1]);
        double vy = (v2[id+2] - v1[id+2]) * (v4[id] - v3[id]) - (v2[id] - v1[id]) * (v4[id+2] - v3[id+2]);
        double vz = (v2[id] - v1[id]) * (v4[id+1] - v3[id+1]) - (v2[id+1] - v1[id+1]) * (v4[id] - v3[id]);

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    /// CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProductmixedID(double *v,
                              double const *v1,
                              double const *v2,
                              double const *v3,
                              double const *v4, int const id1, int const id2,int const id3, int const id4) {
        double vx =
                (v2[id2+1] - v1[id1+1]) * (v4[id4+2] - v3[id3+2]) - (v2[id2+2] - v1[id1+2]) * (v4[id4+1] - v3[id3+1]);
        double vy = (v2[id2+2] - v1[id1+2]) * (v4[id4] - v3[id3]) - (v2[id2] - v1[id1]) * (v4[id4+2] - v3[id3+2]);
        double vz = (v2[id2] - v1[id1]) * (v4[id4+1] - v3[id3+1]) - (v2[id2+1] - v1[id1+1]) * (v4[id4] - v3[id3]);

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3), but with v -> v+d*p. Returns a 3d vector.
    /// ARRAY VERSION
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
    //ARRAY VERSION
    inline void vectorProductStretched(double *v,
                                       double const *v1,
                                       double const *p1,
                                       double const *v2,
                                       double const *p2,
                                       double const *v3,
                                       double const *p3,
                                       double const *v4,
                                       double const *p4,
                                       double d) {
        double vx =
                ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]))
                - ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]));

        double vy =
                ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]))
                - ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]));

        double vz =
                ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]))
                - ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]));
        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    }
    ///CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProductStretched(double *v,
                                       double const *v1,
                                       double const *p1,
                                       double const *v2,
                                       double const *p2,
                                       double const *v3,
                                       double const *p3,
                                       double const *v4,
                                       double const *p4,
                                       double d, int const id) {
        double vx =
                ((v2[id+1] + d * p2[id+1]) - (v1[id+1] + d * p1[id+1])) *
                ((v4[id+2] + d * p4[id+2]) - (v3[id+2] + d * (p3[id+2])))
                - ((v2[id+2] + d * (p2[id+2])) - (v1[id+2] + d * (p1[id+2]))) *
                  ((v4[id+1] + d * (p4[id+1])) - (v3[id+1] + d * (p3[id+1])));

        double vy =
                ((v2[id+2] + d * p2[id+2]) - (v1[id+2] + d * p1[id+2])) *
                (v4[id] + d * p4[id] - (v3[id] + d * p3[id]))
                - ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                  ((v4[id+2] + d * (p4[id+2])) - (v3[id+2] + d * p3[id+2]));

        double vz =
                ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                ((v4[id+1] + d * p4[id+1]) - (v3[id+1] + d * p3[id+1]))
                - ((v2[id+1] + d * p2[id+1]) - (v1[id+1] + d * p1[id+1])) *
                  ((v4[id] + d * p4[id]) - (v3[id] + d * p3[id]));

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    ///CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProductStretchedmixedID(double *v,
                                       double const *v1,
                                       double const *p1,
                                       double const *v2,
                                       double const *p2,
                                       double const *v3,
                                       double const *p3,
                                       double const *v4,
                                       double const *p4,
                                       double d, int const id1, int const id2, int const id3, int const id4) {
        double vx =
                ((v2[id2+1] + d * p2[id2+1]) - (v1[id1+1] + d * p1[id1+1])) *
                ((v4[id4+2] + d * p4[id4+2]) - (v3[id3+2] + d * (p3[id3+2])))
                - ((v2[id2+2] + d * (p2[id2+2])) - (v1[id1+2] + d * (p1[id1+2]))) *
                  ((v4[id4+1] + d * (p4[id4+1])) - (v3[id3+1] + d * (p3[id3+1])));

        double vy =
                ((v2[id2+2] + d * p2[id2+2]) - (v1[id1+2] + d * p1[id1+2])) *
                (v4[id4] + d * p4[id4] - (v3[id3] + d * p3[id3]))
                - ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                  ((v4[id4+2] + d * (p4[id4+2])) - (v3[id3+2] + d * p3[id3+2]));

        double vz =
                ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                ((v4[id4+1] + d * p4[id4+1]) - (v3[id3+1] + d * p3[id3+1]))
                - ((v2[id2+1] + d * p2[id2+1]) - (v1[id1+1] + d * p1[id1+1])) *
                  ((v4[id4] + d * p4[id4]) - (v3[id3] + d * p3[id3]));

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z]. Returns a 3d vector.

    inline vector<double> crossProduct(const vector<double> &v1,
                                       const vector<double> &v2) {
        vector<double> v;

        double vx = v1[1] * v2[2] - v1[2] * v2[1];
        double vy = v1[2] * v2[0] - v1[0] * v2[2];
        double vz = v1[0] * v2[1] - v1[1] * v2[0];

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;
    };
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z].
    /// ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void crossProduct(double *cp,
                             double const *v1,
                             double const *v2) {
        cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
        cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
        cp[2] = v1[0] * v2[1] - v1[1] * v2[0];
    };

    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z]. Returns a 3d vector.
    inline vector<double> crossProductStretched(const vector<double> &v1,
                                                const vector<double> &p1,
                                                const vector<double> &v2,
                                                const vector<double> &p2,
                                                double d) {
        vector<double> v;

        double vx = (v1[1] + d * p1[1]) * (v2[2] + d * p2[2]) - (v1[2] + d * p1[2]) * (v2[1] + d * p2[1]);
        double vy = (v1[2] + d * p1[2]) * (v2[0] + d * p2[0]) - (v1[0] + d * p1[0]) * (v2[2] + d * p2[2]);
        double vz = (v1[0] + d * p1[0]) * (v2[1] + d * p2[1]) - (v1[1] + d * p1[1]) * (v2[0] + d * p2[0]);

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;
    };

    /// Projection of a new point based on a given direction and starting point
    inline vector<double> nextPointProjection(const vector<double> &coordinate,
                                              double d, const vector<double> &tau) {
        vector<double> v;
        v.push_back(coordinate[0] + d * tau[0]);
        v.push_back(coordinate[1] + d * tau[1]);
        v.push_back(coordinate[2] + d * tau[2]);
        return v;
    }

    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v1| = alpha.
    inline vector<double> midPointCoordinate(const vector<double> &v1,
                                             const vector<double> &v2, double alpha) {
        vector<double> v;
        v.push_back(v1[0] * (1.0 - alpha) + alpha * v2[0]);
        v.push_back(v1[1] * (1.0 - alpha) + alpha * v2[1]);
        v.push_back(v1[2] * (1.0 - alpha) + alpha * v2[2]);
        return v;
    }
    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v1| = alpha. ARRAY VERSION

    inline void midPointCoordinate(double *v, double const *v1, double const *v2, double alpha) {

        *(v) = ((*v1) * (1.0 - alpha) + alpha * (*(v2)));
        *(v + 1) = (*(v1 + 1) * (1.0 - alpha) + alpha * (*(v2 + 1)));
        *(v + 2) = (*(v1 + 2) * (1.0 - alpha) + alpha * (*(v2 + 2)));
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void midPointCoordinate(double *v, double const *v1, double const *v2, double alpha, int id) {
        v[0] = v1[id] * (1.0 - alpha) + alpha * v2[id];
        v[1] = v1[id + 1] * (1.0 - alpha) + alpha * v2[id + 1];
        v[2] = v1[id + 2] * (1.0 - alpha) + alpha * v2[id + 2];
    }

    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v| = alpha, but with x-d*p coordinates
    inline vector<double> midPointCoordinateStretched(const vector<double> &v1,
                                                      const vector<double> &p1,
                                                      const vector<double> &v2,
                                                      const vector<double> &p2,
                                                      double alpha, double d) {

        vector<double> v;
        v.push_back((v1[0] + d * p1[0]) * (1.0 - alpha) + alpha * (v2[0] + d * p2[0]));
        v.push_back((v1[1] + d * p1[1]) * (1.0 - alpha) + alpha * (v2[1] + d * p2[1]));
        v.push_back((v1[2] + d * p1[2]) * (1.0 - alpha) + alpha * (v2[2] + d * p2[2]));
        return v;
    }

    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v| = alpha, but with x-d*p coordinates
    /// ARRAY VERSION
    inline void midPointCoordinateStretched(double *v,
                                            double const *v1,
                                            double const *p1,
                                            double const *v2,
                                            double const *p2,
                                            double alpha, double d) {

        v[0] = (v1[0] + d * p1[0]) * (1.0 - alpha) + alpha * (v2[0] + d * p2[0]);
        v[1] = ((v1[1] + d * p1[1]) * (1.0 - alpha) + alpha * (v2[1] + d * p2[1]));
        v[2] = ((v1[2] + d * p1[2]) * (1.0 - alpha) + alpha * (v2[2] + d * p2[2]));
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void midPointCoordinateStretched(double *v,
                                            double *v1,
                                            double *p1,
                                            double *v2,
                                            double *p2,
                                            double alpha, double d, int id) {

        v[0] = (v1[id] + d * p1[id]) * (1.0 - alpha) + alpha * (v2[id] + d * p2[id]);
        v[1] = ((v1[id + 1] + d * p1[id + 1]) * (1.0 - alpha) + alpha * (v2[id + 1] + d * p2[id + 1]));
        v[2] = ((v1[id + 2] + d * p1[id + 2]) * (1.0 - alpha) + alpha * (v2[id + 2] + d * p2[id + 2]));
//        printf("%f \n",(v1[id] + d * p1[id]) * (1.0 - alpha) + alpha * (v2[id] + d * p2[id]));
    }

    /// Returns true if two vectors (p1->p2 and p3->p4) are parallel
    inline bool areParallel(const vector<double> &p1, const vector<double> &p2,
                            const vector<double> &p3, const vector<double> &p4) {

        auto v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        auto v2 = {p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]};

        return areEqual(magnitude(crossProduct(v1, v2)), 0.0);
    }

    /// ARRAY VERSION
    inline bool areParallel(double const *p1, double const *p2,
                            double const *p3, double const *p4) {

        double *v1 = new double[3];
        double *v2 = new double[3];
        double *cp = new double[3];


        v1[0] = p2[0] - p1[0];
        v1[1] = p2[1] - p1[1];
        v1[2] = p2[2] - p1[2];

        v2[0] = p4[0] - p3[0];
        v2[1] = p4[1] - p3[1];
        v2[2] = p4[2] - p3[2];

        crossProduct(cp, v1, v2);

        auto retVal = areEqual(magnitude(cp), 0.0);
        delete [] v1;
        delete [] v2;
        delete [] cp;

        return retVal;
    }

    /// CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline bool areParallel(double const *p1, double const *p2,
                            double const *p3, double const *p4, int id) {

        double v1[3];
        double v2[3];
        double cp[3];


        v1[0] = p2[id] - p1[id];
        v1[1] = p2[id + 1] - p1[id + 1];
        v1[2] = p2[id + 2] - p1[id + 2];

        v2[0] = p4[id] - p3[id];
        v2[1] = p4[id + 1] - p3[id + 1];
        v2[2] = p4[id + 2] - p3[id + 2];

        crossProduct(cp, v1, v2);
//        printf("cp %d %f %f %f\n", id, cp[0], cp[1], cp[2]);
        auto retVal = areEqual(magnitude(cp), 0.0);
//        printf("aE %d %d \n",id, retVal);
//        delete [] v1;
//        delete [] v2;
//        delete [] cp;

        return retVal;
    }

    /// Returns true if two vectors (p1->p2 and p3->p4) are in the same plane
    /// ARRAY VERSION
    inline bool areInPlane(double const *p1, double const *p2,
                           double const *p3, double const *p4) {

        double *v1 = new double[3];
        double *v2 = new double[3];
        double *v3 = new double[3];
        double *cp = new double[3];

        *(v1) = *(p2) - *(p1);
        *(v1 + 1) = *(p2 + 1) - *(p1 + 1);
        *(v1 + 2) = *(p2 + 2) - *(p1 + 2);

        *(v2) = *(p3) - *(p1);
        *(v2 + 1) = *(p3 + 1) - *(p1 + 1);
        *(v2 + 2) = *(p3 + 2) - *(p1 + 2);

        *(v3) = *(p4) - *(p1);
        *(v3 + 1) = *(p4 + 1) - *(p1 + 1);
        *(v3 + 2) = *(p4 + 2) - *(p1 + 2);

        crossProduct(cp, v1, v2);
//        std::cout<<*(p2)<<endl;
        auto retVal = areEqual(dotProduct(v3, cp), 0.0);
        delete [] v1;
        delete [] v2;
        delete [] v3;
        delete [] cp;

        return retVal;
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline bool areInPlane(double const *p1, double const *p2,
                           double const *p3, double const *p4, int const id) {

        double v1[3];
        double v2[3];
        double v3[3];
        double cp[3];

        v1[0] = p2[id] - p1[id];
        v1[1] = p2[id + 1] - p1[id + 1];
        v1[2] = p2[id + 2] - p1[id + 2];

        v2[0] = p3[id] - p1[id];
        v2[1] = p3[id + 1] - p1[id + 1];
        v2[2] = p3[id + 2] - p1[id + 2];

        v3[0] = p4[id] - p1[id];
        v3[1] = p4[id + 1] - p1[id + 1];
        v3[2] = p4[id + 2] - p1[id + 2];

//        printf("%f %f %f %f %f %f %f %f %f\n",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],v3[0],v3[1],v3[2]);
//        cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
//        cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
//        cp[2] = v1[0] * v2[1] - v1[1] * v2[0];

        crossProduct(cp, v1, v2);
//        auto xxx=dotProduct(v3, cp);
        auto retVal = areEqual(dotProduct(v3, cp), 0.0);
//        delete v1, v2, cp;
        return retVal;
    }

    /// Function to move bead out of plane by specified amount
    inline vector<double> movePointOutOfPlane(const vector<double> &p1,
                                              const vector<double> &p2,
                                              const vector<double> &p3,
                                              const vector<double> &p4,
                                              int i, double d) {
        vector<double> norm;
        vector<double> v1;
        vector<double> v2;

        //get plane
        v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        v2 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};

        norm = normalizeVector(crossProduct(v1, v2));

        //move bead 1
        if (i == 1) {
            vector<double> newP1;
            newP1.push_back(p1[0] + norm[0] * d);
            newP1.push_back(p1[1] + norm[1] * d);
            newP1.push_back(p1[2] + norm[2] * d);
            return newP1;
        }

            //move bead 2
        else if (i == 2) {
            vector<double> newP2;
            newP2.push_back(p2[0] + norm[0] * d);
            newP2.push_back(p2[1] + norm[1] * d);
            newP2.push_back(p2[2] + norm[2] * d);
            return newP2;
        }

            //move bead 3
        else if (i == 3) {
            vector<double> newP3;
            newP3.push_back(p3[0] + norm[0] * d);
            newP3.push_back(p3[1] + norm[1] * d);
            newP3.push_back(p3[2] + norm[2] * d);
            return newP3;
        }

            //move bead 4
        else {
            vector<double> newP4;
            newP4.push_back(p4[0] + norm[0] * d);
            newP4.push_back(p4[1] + norm[1] * d);
            newP4.push_back(p4[2] + norm[2] * d);
            return newP4;
        }
    }

    ///IN-PLACE ARRAY VERSION
    inline void movePointOutOfPlane(double *p1,
                                    double *p2,
                                    double *p3,
                                    double *p4,double *newp,
                                    int i, double d) {

        double *norm = new double[3];
        double *v1 = new double[3];
        double *v2 = new double[3];

        //get plane
        v1[0] = p2[0] - p1[0];
        v1[1] = p2[1] - p1[1];
        v1[2] = p2[2] - p1[2];

        v2[0] = p3[0] - p1[0];
        v2[1] = p3[1] - p1[1];
        v2[2] = p3[2] - p1[2];

        crossProduct(norm, v1, v2);
        normalizeVector(norm);

        //move bead 1
        if (i == 1) {

            newp[0] = (p1[0] + norm[0] * d);
            newp[1] = (p1[1] + norm[1] * d);
            newp[2] = (p1[2] + norm[2] * d);
        }

            //move bead 2
        else if (i == 2) {
            newp[0] = (p2[0] + norm[0] * d);
            newp[1] = (p2[1] + norm[1] * d);
            newp[2] = (p2[2] + norm[2] * d);
        }

            //move bead 3
        else if (i == 3) {
            newp[0] = (p3[0] + norm[0] * d);
            newp[1] = (p3[1] + norm[1] * d);
            newp[2] = (p3[2] + norm[2] * d);
        }

            //move bead 4
        else {
            newp[0] = (p4[0] + norm[0] * d);
            newp[1] = (p4[1] + norm[1] * d);
            newp[2] = (p4[2] + norm[2] * d);
        }
        delete [] norm;
        delete [] v1;
        delete [] v2;
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void movePointOutOfPlane(double *p1,
                                    double *p2,
                                    double *p3,
                                    double *p4,
                                    int i, double d, int id) {

        double norm[3];
        double v1[3];
        double v2[3];

        //get plane
        v1[0] = p2[id] - p1[id];
        v1[1] = p2[id + 1] - p1[id + 1];
        v1[2] = p2[id + 2] - p1[id + 2];

        v2[0] = p3[id] - p1[id];
        v2[1] = p3[id + 1] - p1[id + 1];
        v2[2] = p3[id + 2] - p1[id + 2];

        crossProduct(norm, v1, v2);
        normalizeVector(norm);

        //move bead 1
        if (i == 1) {

            p1[id] = (p1[id] + norm[0] * d);
            p1[id + 1] = (p1[id + 1] + norm[1] * d);
            p1[id + 2] = (p1[id + 2] + norm[2] * d);
        }

            //move bead 2
        else if (i == 2) {
            p2[id] = (p2[id] + norm[id] * d);
            p2[id + 1] = (p2[id + 1] + norm[id + 1] * d);
            p2[id + 2] = (p2[id + 2] + norm[id + 2] * d);
        }

            //move bead 3
        else if (i == 3) {
            p3[id] = (p3[0] + norm[0] * d);
            p3[id + 1] = (p3[id + 1] + norm[1] * d);
            p3[id + 2] = (p3[id + 2] + norm[2] * d);
        }

            //move bead 4
        else {
            p4[id] = (p4[id] + norm[0] * d);
            p4[id + 1] = (p4[id + 1] + norm[1] * d);
            p4[id + 2] = (p4[id + 2] + norm[2] * d);
        }
//        delete [] norm;
//        delete [] v1;
//        delete [] v2;
    }

    /// Function to create a initial branching point and direction, given an
    /// initial normal vector and point.
    /// @param l - the distance of the branch from the original point
    /// @param m - the size of the branch projection
    /// @param theta - the angle of branching
    /// @return a vector describing the initial branching direction and point
    tuple<vector<double>, vector<double>> branchProjection(const vector<double> &n,
                                                           const vector<double> &p,
                                                           double l, double m, double theta);


    /// Returns true if two vectors (p1->p2 and p3->p4) are in the same plane
    inline bool areInPlane(const vector<double> &p1, const vector<double> &p2,
                           const vector<double> &p3, const vector<double> &p4) {

        auto v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        auto v2 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};
        auto v3 = {p4[0] - p1[0], p4[1] - p1[1], p4[2] - p1[2]};

        auto cp = crossProduct(v1, v2);

        return areEqual(dotProduct(v3, cp), 0.0);
    }

    inline size_t blockToSmem(int blockSize){return 12 * blockSize * sizeof(double);}
    inline size_t blockToSmemez(int blockSize){return 24 * blockSize * sizeof(double);}
    inline size_t blockToSmemF(int blockSize){return 6 * blockSize * sizeof(double);}
    inline size_t blockToSmemFB(int blockSize){return 9 * blockSize * sizeof(double);}
    inline size_t blockToSmemFB2(int blockSize){return 18 * blockSize * sizeof(double);}
    inline size_t blockToSmemFB3(int blockSize){return 3 * blockSize * sizeof(double);}


    /// Function to move bead out of plane by specified amount
    vector<double> movePointOutOfPlane(const vector<double> &p1,
                                       const vector<double> &p2,
                                       const vector<double> &p3,
                                       const vector<double> &p4,
                                       int i, double d);


    /// Function to create a initial branching point and direction, given an
    /// initial normal vector and point.
    /// @param l - the distance of the branch from the original point
    /// @param m - the size of the branch projection
    /// @param theta - the angle of branching
    /// @return a vector describing the initial branching direction and point
    tuple<vector<double>, vector<double>> branchProjection(const vector<double> &n,
                                                           const vector<double> &p,
                                                           double l, double m, double theta);
    
    
    float delGGenChem(float delGZero, vector<species_copy_t> reacN, vector<int> reacNu, vector<species_copy_t> prodN, vector<int> prodNu);
    
    float delGGenChemI(float delGZero, vector<species_copy_t> reacN, vector<int> reacNu, vector<species_copy_t> prodN, vector<int> prodNu);
    
    float delGDifChem(species_copy_t reacN ,species_copy_t prodN);
    
    
    float delGPolyChem(float delGzero, species_copy_t reacN, string whichWay);
    
    float delGMyoChem(float nh, float rn);

}
#endif
