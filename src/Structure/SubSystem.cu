
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

#include "SubSystem.h"

#include "BoundaryElement.h"
#include "CompartmentGrid.h"
#include "BindingManager.h"
#include "BindingManagerCUDA.h"
#include "CUDAcommon.h"
#include "MathFunctions.h"
using namespace mathfunc;
void SubSystem::updateBindingManagers() {
#ifdef CUDAACCL
    if(SysParams::Chemistry().numFilaments > 1) {
        cout << "CUDA Binding Manager cannot handle more than one type of filaments." << endl;
        exit(EXIT_FAILURE);
    }
    numpairs_vec.clear();
    possibleBindings_vec.clear();
    gpu_possibleBindings_vec.clear();
//   auto x = CMonomer::_numBSpecies;
//    auto var = SysParams::Chemistry().bmanagerdistances;
    Compartment* C0 = _compartmentGrid->getCompartments()[0];
    for(auto &manager : C0->getFilamentBindingManagers()) {

        LinkerBindingManager *lManager;
        MotorBindingManager *mManager;
        BranchingManager *bManager;
        auto cylcylnlvars = CUDAcommon::getCylCylNLvars();
        auto coord = cylcylnlvars.gpu_coord;
        auto beadSet = cylcylnlvars.gpu_beadSet;
        auto cylID = cylcylnlvars.gpu_cylID;
        auto filType = cylcylnlvars.gpu_filType;
        auto filID = cylcylnlvars.gpu_filID;
        int *cmpID = cylcylnlvars.gpu_cmpID;

        if (gpu_params == NULL) {
            int params[3];
            params[0] = SysParams::Chemistry().numBindingSites[0];//filType dependant
            params[1] = 0;//filType dependant
            params[2] = SysParams::Geometry().cylinderNumMon[0];//filType dependant.
            CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)), "cuda data transfer",
                                    "SubSystem.cu");
            CUDAcommon::handleerror(cudaMemcpy(gpu_params, params, 3 * sizeof(int), cudaMemcpyHostToDevice));
        }
        if (gpu_bindingSites == NULL) {
            auto bindingSites = SysParams::Chemistry().bindingSites[0];//filType dependant
            int cpu_bindingSites[bindingSites.size()];
            int iii = 0;
            for (auto bs:bindingSites)
            {cpu_bindingSites[iii] = int(bs); iii++;}
            CUDAcommon::handleerror(cudaMalloc((void **) &gpu_bindingSites, bindingSites.size() *
                                            sizeof(int)), "cuda data transfer", "SubSystem.cu");
            CUDAcommon::handleerror(cudaMemcpy(gpu_bindingSites, cpu_bindingSites, bindingSites.size() *  sizeof
                                    (int), cudaMemcpyHostToDevice));
        }
        //Linker
        if ((lManager = dynamic_cast<LinkerBindingManager *>(manager.get()))) {
            cudaStream_t  s;
            if(numbindmgrs + 1 > strvec.size() )
            { cudaStreamCreate(&s); strvec.push_back(s);}
            else
                s = strvec.at(numbindmgrs);
            numbindmgrs++;
            int *cmon_state_linker = cylcylnlvars.gpu_cmon_state_linker;
            //2. Assign optimal blocks and threads
            vector<int> blocksnthreads;
            int blockSize;   // The launch configurator returned block size
            int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
            int nint = lManager->getNLsize();
            int *NL = lManager->getNLCUDA();
            int *numNLpairs = lManager->getNLsizeCUDA();
            int *numpairs = lManager->getpossiblebindingssizeCUDA();
            double *params2 = lManager->getdistancesCUDA();
            std::cout<<"Total Linker NL size "<<nint<<endl;
//            int *numpairs, test[1];test[0] = 0;
//            CUDAcommon::handleerror(cudaMalloc((void **) &numpairs, sizeof(int)), "cuda data transfer", "SubSystem.cu");
//            CUDAcommon::handleerror(cudaMemcpy(numpairs, test, sizeof(int), cudaMemcpyHostToDevice));
            if (nint > 0) {
                int *gpu_possibleBindings;

                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_possibleBindings, SysParams::Chemistry()
                        .numBindingSites[0] * 5 * nint * sizeof(int)));
                cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                                   updateAllPossibleBindingsCUDA, 0, 0);
                blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
                blocksnthreads.push_back(blockSize);
                std::cout << "Linker blocks and threads " << blocksnthreads.at(0) << " " << blocksnthreads.at(1)
                          << endl;
                resetintvariableCUDA<<<1,1,0,s>>>(numpairs);
                updateAllPossibleBindingsCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                      (coord, beadSet, cylID, filID, filType, cmpID, NL, numNLpairs, numpairs,
                                                                              gpu_params, params2, gpu_possibleBindings, cmon_state_linker,
                                                                              gpu_bindingSites);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
                //copy back to CPU
                int *cpu_numpairs, *possibleBindings;
                CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                        "Subsystem.cu");
                CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                        s),"Copy", "Subsystem.cu");
                CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
                CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
                CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
                        cudaHostAllocDefault),"Copy", "Subsystem.cu");
                CUDAcommon::handleerror(cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] *
                sizeof(int), cudaMemcpyDeviceToHost, s),"Copy", "Subsystem.cu");
                numpairs_vec.push_back(cpu_numpairs);
                possibleBindings_vec.push_back(possibleBindings);
                gpu_possibleBindings_vec.push_back(gpu_possibleBindings);
//                int cpu_numpairs[1];
//                cudaMemcpy(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost);
//                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
//                int possibleBindings[5 * cpu_numpairs[0]];
//                cudaMemcpy(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
//                           cudaMemcpyDeviceToHost);
            }
        }
        //Motor
        else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
            cudaStream_t  s;
            if(numbindmgrs + 1 > strvec.size() )
            { cudaStreamCreate(&s); strvec.push_back(s);}
            else
                s = strvec.at(numbindmgrs);
            numbindmgrs++;
            int *cmon_state_linker = cylcylnlvars.gpu_cmon_state_motor;
            //2. Assign optimal blocks and threads
            vector<int> blocksnthreads;
            int blockSize;   // The launch configurator returned block size
            int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
            int nint = mManager->getNLsize();
            int *NL = mManager->getNLCUDA();
            int *numNLpairs = mManager->getNLsizeCUDA();
            int *numpairs = mManager->getpossiblebindingssizeCUDA();
            double *params2 = mManager->getdistancesCUDA();
            std::cout<<"Total Motor NL size "<<nint<<endl;
            if (nint > 0) {
                int *gpu_possibleBindings;

                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_possibleBindings, SysParams::Chemistry()
                        .numBindingSites[0] * 5 * nint * sizeof(int)));
                cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                                   updateAllPossibleBindingsCUDA, 0, 0);
                blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
                blocksnthreads.push_back(blockSize);
                std::cout << "Motor blocks and threads " << blocksnthreads.at(0) << " " << blocksnthreads.at(1)
                          << endl;
                resetintvariableCUDA << < 1, 1, 0, s >> > (numpairs);
                updateAllPossibleBindingsCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                      (coord, beadSet, cylID, filID, filType, cmpID, NL, numNLpairs, numpairs,
                                                                              gpu_params, params2, gpu_possibleBindings, cmon_state_linker,
                                                                              gpu_bindingSites);

                //copy back to CPU
                int *cpu_numpairs, *possibleBindings;
                CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                        "Subsystem.cu");
                CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                        s),"Copy", "Subsystem.cu");
                CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
                CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
                CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
                                                      cudaHostAllocDefault),"Copy", "Subsystem.cu");
                CUDAcommon::handleerror(cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] *
                                        sizeof(int), cudaMemcpyDeviceToHost, s),"Copy", "Subsystem.cu");
                numpairs_vec.push_back(cpu_numpairs);
                possibleBindings_vec.push_back(possibleBindings);
                gpu_possibleBindings_vec.push_back(gpu_possibleBindings);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
//                CUDAcommon::handleerror(cudaDeviceSynchronize());
//                //copy back to CPU
//                int cpu_numpairs[1];
//                cudaMemcpy(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost);
//                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
//                int possibleBindings[5 * cpu_numpairs[0]];
//                cudaMemcpy(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
//                           cudaMemcpyDeviceToHost);

            }
            }
        //Brancher
        else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {

            }

    }
    //Synchronize streams
    for(auto s:strvec) CUDAcommon::handleerror(cudaStreamSynchronize(s),"stream sync","SubsSystem.cu");
    //clear all possible bindings.
    for(auto c:_compartmentGrid->getCompartments()){
        for(auto &Mgr : c->getFilamentBindingManagers()) {
            Mgr->clearpossibleBindings();
        }
    }
    //Assign to respective possible bindings.
    int count = 0;
    for(auto &manager : C0->getFilamentBindingManagers()) {
        LinkerBindingManager *lManager;
        MotorBindingManager *mManager;
        BranchingManager *bManager;
        //Linkers
        if ((lManager = dynamic_cast<LinkerBindingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[5* i];
                int cIndex = possibleBindings[5 * i +1];
                short cbs = short(possibleBindings[5 * i + 2]);
                int cnIndex = possibleBindings[5 * i +3];
                short cnbs = short(possibleBindings[5 * i + 4]);
                auto cylinder = Cylinder::getCylinders()[cIndex];
                auto ncylinder = Cylinder::getCylinders()[cnIndex];
                //get the compartment.
                Compartment* cmp = GController::getCompartment(cID);
                //get corresponding binding manager
                for(auto &cmanager : cmp->getFilamentBindingManagers()) {
                    if ((lManager = dynamic_cast<LinkerBindingManager *>(cmanager.get()))) {
                        auto t1 = tuple<CCylinder*, short>(cylinder->getCCylinder(), cbs);
                        auto t2 = tuple<CCylinder*, short>(ncylinder->getCCylinder(), cnbs);
                        cmanager->appendpossibleBindings(t1,t2);
                    }
                }
            }
            count++;
        }
        //MOTORS
        else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[5* i];
                int cIndex = possibleBindings[5 * i +1];
                short cbs = short(possibleBindings[5 * i + 2]);
                int cnIndex = possibleBindings[5 * i +3];
                short cnbs = short(possibleBindings[5 * i + 4]);
                auto cylinder = Cylinder::getCylinders()[cIndex];
                auto ncylinder = Cylinder::getCylinders()[cnIndex];
                //get the compartment
                Compartment* cmp = GController::getCompartment(cID);
                //get corresponding binding manager
                for(auto &cmanager : cmp->getFilamentBindingManagers()) {
                    if ((mManager = dynamic_cast<MotorBindingManager *>(cmanager.get()))) {
                        auto t1 = tuple<CCylinder*, short>(cylinder->getCCylinder(), cbs);
                        auto t2 = tuple<CCylinder*, short>(ncylinder->getCCylinder(), cnbs);
                        cmanager->appendpossibleBindings(t1,t2);
                    }
                }
            }
            count++;
        }
        else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {

        }
    }
    for(auto gpb:gpu_possibleBindings_vec)
        CUDAcommon::handleerror(cudaFree(gpb),"cudaFree","SubSystem.cu");
    for(auto pb:possibleBindings_vec)
        CUDAcommon::handleerror(cudaFreeHost(pb),"cudaFree","SubSystem.cu");
    for(auto np:numpairs_vec)
        CUDAcommon::handleerror(cudaFreeHost(np),"cudaFree","SubSystem.cu");

    //cudaFree
    endresetCUDA();
#else

    for(auto C : _compartmentGrid->getCompartments()) {

        for(auto &manager : C->getFilamentBindingManagers())
            manager->updateAllPossibleBindings();
    }
#endif

//    int l =0; int m = 0;int b = 0;
//    for(auto C: _compartmentGrid->getCompartments()){
//        for(auto &manager : C->getFilamentBindingManagers()) {
//            LinkerBindingManager *lManager;
//            MotorBindingManager *mManager;
//            BranchingManager *bManager;
//            if ((lManager = dynamic_cast<LinkerBindingManager *>(manager.get()))) {
//                l += lManager->numBindingSites();
//            } else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
//                m += mManager->numBindingSites();
//            } else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {
//                b += bManager->numBindingSites();
//            }
//        }
//    }
//    std::cout<<"L M B "<<l<<" "<<m<<" "<<b<<endl;
//    std::cout<<endl;
}

