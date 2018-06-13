//
// Created by aravind on 9/18/17.
//

#ifndef CUDA_VEC_CUDACOMMON_H
#define CUDA_VEC_CUDACOMMON_H

#ifdef CUDAACCL
#include <vector>
#include <list>
#include <src/Mechanics/ForceField/Filament/FilamentStretchingHarmonic.h>
#include <src/Mechanics/ForceField/Filament/FilamentBendingHarmonic.h>
#include <src/Mechanics/ForceField/Filament/FilamentBendingCosine.h>
#include <src/Mechanics/ForceField/Linker/LinkerStretchingHarmonic.h>
#include <src/Mechanics/ForceField/MotorGhost/MotorGhostStretchingHarmonic.h>
#include <src/Mechanics/ForceField/Volume/CylinderExclVolRepulsion.h>
#include <src/Mechanics/ForceField/Branching/BranchingStretchingHarmonic.h>
#include <src/Mechanics/ForceField/Branching/BranchingBendingCosine.h>
#include <src/Mechanics/ForceField/Branching/BranchingDihedralCosine.h>
#include <src/Mechanics/ForceField/Branching/BranchingPositionCosine.h>
#include <src/Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.h>

#include "common.h"
#include "string.h"
#include "MathFunctions.h"
using namespace mathfunc;
struct CUDAvars {
    double * gpu_force = NULL;
    double * gpu_forceAux = NULL;
    double * gpu_forceAuxP = NULL;
    double * gpu_coord =  NULL;
    double * gpu_lambda = NULL;
    float vectorize = 0.0;
    double * gpu_energy = NULL;
    bool * gpu_btstate = NULL;
    vector<cudaStream_t*> streamvec;
    vector<cudaEvent_t> eventvec;
    int* culpritID = NULL;
    int* gculpritID = NULL;
    char* gculpritFF = NULL;
    char* gculpritinteraction = NULL;
    char* culpritFF = NULL;
    char* culpritinteraction = NULL;
    size_t memincuda = 0;
    bool conservestreams = true;
//    cudaEvent_t *event;

//    float Ccforce = 0.0;
//    float Scforce = 0.0;
//    unsigned int  gpu_sharedMem = 0;
//    unsigned int  gpu_globalMem = 0;
//    int * motorparams;
};

struct CylCylNLvars {
    double* gpu_coord;
    double* gpu_coord_com;
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
};

class CUDAcommon{
public:
    static CUDAvars cudavars;
    static CylCylNLvars cylcylnlvars;
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

};
#endif
#endif
//CUDA_VEC_CUDACOMMON_H
