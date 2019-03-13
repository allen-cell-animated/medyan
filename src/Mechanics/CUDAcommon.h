//
// Created by aravind on 9/18/17.
//

#ifndef CUDA_VEC_CUDACOMMON_H
#define CUDA_VEC_CUDACOMMON_H
#include <vector>
#include <list>
#include <FilamentStretchingHarmonic.h>
#include <FilamentBendingHarmonic.h>
#include <FilamentBendingCosine.h>
#include <LinkerStretchingHarmonic.h>
#include <MotorGhostStretchingHarmonic.h>
#include <CylinderExclVolRepulsion.h>
#include <BranchingStretchingHarmonic.h>
#include <BranchingBendingCosine.h>
#include <BranchingDihedralCosine.h>
#include <BranchingPositionCosine.h>
#include <BoundaryCylinderRepulsionExp.h>
#include "CCylinder.h"
#include "common.h"
#include "string.h"
#include "MathFunctions.h"
#include "dist_driver.h"
using namespace mathfunc;
struct bin{
    int binID;
    floatingpoint bincoord[3];
    int neighbors[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                          -1,-1,-1,-1,-1,-1};
    int binstencilID[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                          -1,-1,-1,-1,-1,-1};
};
    struct cylinder{
        int filamentID = -1;
        int filamentposition = -1;
        int bindices[2];
        int cmpID = -1;
        long cindex = -1;
        floatingpoint coord[3];
        short type = -1;
        int ID = -1;
        int availbscount = -1;
    };
struct SERLvars{
    floatingpoint *coord = NULL;
    cylinder *cylindervec = NULL;
    CCylinder **ccylindervec = NULL;
    Cylinder **cylinderpointervec = NULL;
    uint N = 6000;

};
template <uint D>
dist::dOut<D> SIMDoutvar(const uint dim, uint N1, std::initializer_list<float> params) {

    if (dim == 1) {
        dist::dOut<1> out_serialdim(N1, params);
        return out_serialdim;
    }
    else if (dim == 2) {
        dist::dOut<2> out_serialdim(N1, params);
        return out_serialdim;
    }
    else if (dim == 3){
        dist::dOut<3> out_serialdim(N1, params);
        return out_serialdim;
    }
}
#if defined(CUDAACCL) || defined(CUDATIMETRACK)
struct CUDAvars {
    floatingpoint * gpu_force = NULL;
    floatingpoint * gpu_forceAux = NULL;
    floatingpoint * gpu_forceAuxP = NULL;
    floatingpoint * gpu_coord =  NULL;
    floatingpoint * gpu_lambda = NULL;
    float vectorize = 0.0;
    floatingpoint * gpu_energy = NULL;
    bool * gpu_btstate = NULL;
    cylinder* gpu_cylindervec = NULL;
#ifdef CUDAACCL
    vector<cudaStream_t*> streamvec;
    vector<cudaEvent_t> eventvec;
#endif
    int* culpritID = NULL;
    int* gculpritID = NULL;
    char* gculpritFF = NULL;
    char* gculpritinteraction = NULL;
    char* culpritFF = NULL;
    char* culpritinteraction = NULL;
    size_t memincuda = 0;
    bool conservestreams = true;
    int offset_E = 0;
    floatingpoint *gpu_energyvec = NULL;
    vector<bool*> backtrackbools;
//    cudaEvent_t *event;

//    float Ccforce = 0.0;
//    float Scforce = 0.0;
//    unsigned int  gpu_sharedMem = 0;
//    unsigned int  gpu_globalMem = 0;
//    int * motorparams;
};

struct CylCylNLvars {
    floatingpoint* gpu_coord;
    floatingpoint* gpu_coord_com;
    int * gpu_beadSet;
    int *gpu_cylID;
    int *gpu_filID;
    int *gpu_filType;
    int *gpu_cmpID;
    int *gpu_fvecpos;
//    int *gpu_cylstate;
    int *gpu_cmon_state_brancher;
    int *gpu_cmon_state_linker;
    int *gpu_cmon_state_motor;
//    int *gpu_cylvecpospercmp;

    bin* bins;

};

struct SERLtime {
    floatingpoint TvectorizeFF = 0.0;
    floatingpoint TcomputeE = 0.0;
    floatingpoint TcomputeEiter = 0.0;
    int Ecount = 0;
    floatingpoint TcomputeF= 0.0;
    floatingpoint Tlambda = 0.0;
    floatingpoint TshiftGrad= 0.0;
    floatingpoint TmaxF= 0.0;
    floatingpoint Tcalculatedot = 0.0;
    vector<floatingpoint>TvecvectorizeFF;
    vector<floatingpoint>TveccomputeE;
    vector<floatingpoint>TveccomputeF;
    vector<floatingpoint>Tlambdavec;
    vector<floatingpoint>Tlambdap;
    vector<floatingpoint>Tlambdapcount;
};

struct CUDAtime {
    floatingpoint TvectorizeFF = 0.0;
    floatingpoint TcomputeE = 0.0;
    floatingpoint TcomputeEiter = 0.0;
    int Ecount = 0;
    floatingpoint TcomputeF= 0.0;
    floatingpoint Tlambda = 0.0;
    floatingpoint TshiftGrad= 0.0;
    floatingpoint TmaxF= 0.0;
    floatingpoint Tcalculatedot = 0.0;
    floatingpoint Tstartmin = 0.0;
    vector<floatingpoint>TvecvectorizeFF;
    vector<floatingpoint>TveccomputeE;
    vector<floatingpoint>TveccomputeF;
    vector<floatingpoint>Tlambdavec;
    vector<floatingpoint>Tlambdap;
    vector<floatingpoint>Tlambdapcount;
};
#endif
class CUDAcommon{
public:
    static SERLvars serlvars;
    static const SERLvars& getSERLvars(){return serlvars;}
#ifdef CUDAACCL
    static CUDAvars cudavars;
    static CylCylNLvars cylcylnlvars;
    static SERLtime serltime;
    static CUDAtime cudatime;
    static const CUDAvars& getCUDAvars(){return cudavars;}
    static const CylCylNLvars& getCylCylNLvars(){return cylcylnlvars;}
    static void handleerror(cudaError_t a){
        if(a !=cudaSuccess){
            cout<<cudaGetErrorString(a)<<endl;
            exit(EXIT_FAILURE);
        }
    }

    static void handleerror(cudaError_t a, std::string tag1, std::string tag2){
//        cout<<CUDAcommon::getCUDAvars().culpritFF[0]<<CUDAcommon::getCUDAvars().culpritFF[1]<<endl;
        if(a == cudaErrorAssert){
            cudaDeviceSynchronize();
            auto culpritFF = CUDAcommon::getCUDAvars().culpritFF;
            auto culpritinteraction = CUDAcommon::getCUDAvars().culpritinteraction;

            if(strcmp(culpritinteraction, "Filament Stretching Harmonic")==0)
                FilamentStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Filament Bending Harmonic")==0)
                FilamentBendingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Filament Bending Cosine")==0)
                FilamentBendingCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Linker Stretching Harmonic")==0)
                LinkerStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Motor Stretching Harmonic")==0)
                MotorGhostStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Cylinder Excluded Volume")==0)
                CylinderExclVolRepulsion::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Stretching Harmonic")==0)
                BranchingStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Bending Cosine")==0)
                BranchingBendingCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Dihedral Cosine")==0)
                BranchingDihedralCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Position Cosine")==0)
                BranchingPositionCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Boundary Cylinder Repulsion Exp")==0)
                BoundaryCylinderRepulsionExp::checkforculprit();
            else{
                cout<<"unknown assert error. Check code. Exiting.."<<endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(a !=cudaSuccess && a != cudaErrorAssert){
            cout<<cudaGetErrorString(a)<<endl;
            cout<<"ERROR! "<<tag1<<"error in "<<tag2<<". Check vectors. Exiting.."<<endl;
            exit(EXIT_FAILURE);
        }
    }

    static void printculprit(std::string tag1, std::string tag2){
        cout << "Energy of system became infinite. Try adjusting minimization parameters." << endl;
        cout << "The culprit was ... " << tag1 << endl;
        cout << "Culprit interaction = " << tag2 << endl;
    }
#endif
};
#endif
//CUDA_VEC_CUDACOMMON_H
