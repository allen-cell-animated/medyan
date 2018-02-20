
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

#include <cmath>
#include <vector>
#include <math.h>

#include "MathFunctions.h"
#include "Rand.h"

namespace mathfunc {

    tuple<vector<double>, vector<double>> branchProjection(const vector<double>& n,
                                                           const vector<double>& p,
                                                           double l, double m, double theta){
        //get random permutation from p
        vector<double> r = {p[0] + Rand::randDouble(-1, 1),
                            p[1] + Rand::randDouble(-1, 1),
                            p[2] + Rand::randDouble(-1, 1)};

        //construct vector z which is r-p
        auto z = twoPointDirection(p, r);

        //construct u and v, which creates an orthogonal set n, u, v
        auto u = crossProduct(n, z);
        auto v = crossProduct(n, u);

        normalize(u); normalize(v);

        //find random point on circle defining the branching point
        double thetaRandom = Rand::randDouble(0, 2*M_PI);
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
        bp2.push_back(newP[0] + newL * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp2.push_back(newP[1] + newL * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp2.push_back(newP[2] + newL * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));

        //get direction
        auto direction = twoPointDirection(bp1, bp2);

        return tuple<vector<double>, vector<double>>(direction, bp1);
    }
     __global__ void addvector(double *U, int *params, double *U_sum, double *U_tot){
        U_sum[0] = 0.0;
        double sum = 0.0;
        for(auto i=0;i<params[1];i++){
            if(U[i] == -1.0 && sum != -1.0){
                U_sum[0] = -1.0;
                U_tot[0] = -1.0;
                sum = -1.0;
                break;
            }
            else
                sum  += U[i];
        }
        U_sum[0] = sum;
        atomicAdd(&U_tot[0], sum);
//        printf("%f %f\n", U_tot[0], sum);

    }
    __global__ void addvector(double *U, int *params, double *U_sum, double *U_tot, int *culpritID, char* culpritFF,
                              char* culpritinteraction, char* FF, char* interaction){
        U_sum[0] = 0.0;
        double sum = 0.0;
        for(auto i=0;i<params[1];i++){
            if(U[i] == -1.0 && sum != -1.0){
                U_sum[0] = -1.0;
                U_tot[0] = -1.0;
                sum = -1.0;
            }
            else if(sum != -1.0)
                sum  += U[i];
        }
        U_sum[0] = sum;
//        printf("sum %f\n",sum);
        if(sum != -1.0)
            atomicAdd(&U_tot[0], sum);
        else{
            assert(0);
        }
    }
}