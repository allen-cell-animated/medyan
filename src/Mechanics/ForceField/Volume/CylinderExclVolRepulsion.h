
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

#ifndef MEDYAN_CylinderExclVolRepulsion_h
#define MEDYAN_CylinderExclVolRepulsion_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// CylinderExclVolume template.
class CylinderExclVolRepulsion {
    
public:
    double energy(double *coord, double *f, int *beadSet, double *krep);
    
    double energy(double *coord, double *f, int *beadSet, double *krep, double d);
    
    void forces(double *coord, double *f, int *beadSet, double *krep);

#ifdef CUDAACCL
    double energy(double *coord, double *f, int *beadSet, double *krep, int *params, vector<int> blocksnthreads);
#endif
#ifdef CROSSCHECK
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls, double d);
    
    void forces(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double Krepuls);
#endif
};
#ifdef CUDAACCL
//__global__
//void saxpy(int n, float a, float *x, float *y)
//{
//    int i = blockIdx.x*blockDim.x + threadIdx.x;
//    if (i < n) y[i] = a*x[i] + y[i];
//}

//__global__ void CUDAenergy(double *coord, double *force, int *beadSet, double *krep) {
////memory needed: 34*THREADSPERBLOCK*sizeof(double)+2*THREADSPERBLOCK*sizeof(int);
//    extern __shared__ double s[];
//    double *coord_image=s;
//    double *c1 = &s[blockDim.x];
//    double *c2 = &s[blockDim.x];
//    double *c3 = &s[blockDim.x];
//    double *c4 = &s[blockDim.x];
//    double *U_i = &s[blockDim.x];
//    double d, invDSquare;
//    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
//    double ATG1, ATG2, ATG3, ATG4;
////    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
////    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;
//    int n = 4;
//    int nint = 10;
//
//    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    if(thread_idx<nint) {
//        for(auto i=0;i<3;i++){
//            c1[i] = coord[3 * beadSet[n * i]];
//            c2[i] = coord[3 * beadSet[n * i + 1]];
//            c3[i] = coord[3 * beadSet[n * i + 2]];
//            c4[i] = coord[3 * beadSet[n * i + 3]];
//        }}
//
//    __syncthreads();
//
////    if(thread_idx<nint) {
////        //check if parallel
////        if(areParallel(c1, c2, c3, c4)) {
////
////            d = twoPointDistance(c1, c3);
////            invDSquare =  1 / (d * d);
////            U_i = krep[i] * invDSquare * invDSquare;
//////            std::cout<<U_i<<endl;
////            if(fabs(U_i) == numeric_limits<double>::infinity()
////               || U_i != U_i || U_i < -1.0) {
////
////                //set culprit and return TODO
////                return -1;
////            }
////            continue;
////        }
////
////        //check if in same plane
////        if(areInPlane(c1, c2, c3, c4)) {
////
////            //slightly move point
////            movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01);
////        }
////
////        a = scalarProduct(c1, c2, c1, c2);
////        b = scalarProduct(c3, c4, c3, c4);
////        c = scalarProduct(c3, c1, c3, c1);
////        d = scalarProduct(c1, c2, c3, c4);
////        e = scalarProduct(c1, c2, c3, c1);
////        F = scalarProduct(c3, c4, c3, c1);
////
////        AA = sqrt(a*c - e*e);
////        BB = sqrt(b*c - F*F);
////
////        CC = d*e - a*F;
////        DD = b*e - d*F;
////
////        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
////        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );
////
////        GG = d*d - a*b - CC;
////        HH = CC + GG - DD;
////        JJ = c*(GG + CC) + e*DD - F*CC;
////
////
////        ATG1 = atan( (a + e)/AA) - atan(e/AA);
////        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
////        ATG3 = atan((F)/BB) - atan((F - b)/BB);
////        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);
////
////        U_i[thread_idx] = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
////    }
////
////    __syncthreads();
//}
#endif

#endif
