
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ForceFieldManager.h"
#include "ForceFieldManagerCUDA.h"

#include "CGMethod.h"
#include "cross_check.h"
#include <algorithm>

#include "Structure/Bead.h"
#include "Structure/Cylinder.h"







ForceField* ForceFieldManager::_culpritForceField = nullptr;

void ForceFieldManager::vectorizeAllForceFields() {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    CUDAcommon::cudatime.TvectorizeFF = 0.0;
    CUDAcommon::cudatime.TvecvectorizeFF.clear();
#endif
#ifdef CUDAACCL
    // PT1 Generate single vector of energies from all FF and add them together.
    //@{
    CUDAcommon::cudavars.offset_E=0.0;
    //@}
#endif

    for (auto &ff : _forceFields)
        ff->vectorize();

#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //reset offset
    if (streamF == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&streamF));
    int nint[1];
    nint[0] = CGMethod::N / 3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_nint, nint, sizeof(int),
                                        cudaMemcpyHostToDevice, streamF));
    CUDAcommon::handleerror(cudaMalloc((void **) &(CUDAcommon::cudavars.gpu_energyvec),
                                       CUDAcommon::cudavars.offset_E * sizeof(floatingpoint)));
    int THREADSPERBLOCK;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    THREADSPERBLOCK = prop.maxThreadsPerBlock;

    blocksnthreads.push_back(CGMethod::N / (3 * THREADSPERBLOCK) + 1);
    if (blocksnthreads[0] == 1) blocksnthreads.push_back(CGMethod::N / 3);
    else blocksnthreads.push_back(THREADSPERBLOCK);

    // PT2 Generate single vector of energies from all FF and add them together.
    //@{
//    std::cout<<"CUDA energy total nint "<<CUDAcommon::cudavars.offset_E<<endl;
    bntaddvec2.clear();
    bntaddvec2 = getaddred2bnt(CUDAcommon::cudavars.offset_E);
    CUDAcommon::handleerror(cudaMalloc((void **) &(CUDAcommon::cudavars.gpu_energyvec), bntaddvec2.at
            (0)*sizeof (floatingpoint)));
    vector<floatingpoint> zerovec(bntaddvec2.at(0));
    fill(zerovec.begin(),zerovec.begin()+bntaddvec2.at(0),0.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, zerovec.data(),
                            bntaddvec2.at(0) * sizeof(floatingpoint), cudaMemcpyHostToDevice,streamF));
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 0,
                                            bntaddvec2.at(0) * sizeof(floatingpoint), streamF));*/

    params.clear();
    params.push_back(CUDAcommon::cudavars.offset_E);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), sizeof(int),
                                       cudaMemcpyHostToDevice, streamF));
    //@}
#endif
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
    std::cout<<"Time total vectorizeFF (s) "<<CUDAcommon::cudatime.TvectorizeFF<<endl;
    std::cout<<"Time split vectorizeFF (s) ";
    for(auto x:CUDAcommon::cudatime.TvecvectorizeFF)
        std::cout<<x<<" ";
    std::cout<<endl;
#endif
}

void ForceFieldManager::cleanupAllForceFields() {

    for (auto &ff : _forceFields)
        ff->cleanup();
#ifdef CUDAACCL
    if (!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(streamF));
    //cleanup energy vector
    // single vector of energies from all FF and add them together.
    //@{
    CUDAcommon::handleerror(cudaFree(CUDAcommon::cudavars.gpu_energyvec));
    CUDAcommon::handleerror(cudaFree(gpu_params));
    //@}
    if (CGMethod::N / 3 > 0) {
        CUDAcommon::handleerror(cudaFree(gpu_nint));
        //Memory alloted
        //@{
//        size_t allocmem = 0;
//        allocmem += sizeof(floatingpoint);
//        auto c = CUDAcommon::getCUDAvars();
//        c.memincuda -= allocmem;
//        CUDAcommon::cudavars = c;
//        std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
//        std::cout<<"Memory allocated 0 . Memory freed "<<allocmem/1024<<endl;
        //@}
        blocksnthreads.clear();
    }
#endif
}

template< bool stretched >
floatingpoint ForceFieldManager::computeEnergy(floatingpoint *coord, bool verbose) const {
    chrono::high_resolution_clock::time_point tbegin, tend;
#ifdef CUDATIMETRACK
//    CUDAcommon::cudatime.TcomputeE = 0.0;
    CUDAcommon::cudatime.TveccomputeE.clear();
    CUDAcommon::cudatime.Ecount++;
//    CUDAcommon::serltime.TcomputeE = 0.0;
    CUDAcommon::serltime.TveccomputeE.clear();
    CUDAcommon::serltime.Ecount++;
#endif
    floatingpoint energy = 0.0;
#ifdef CUDAACCL
#ifdef CUDA_INDIVIDUAL_ESUM
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (floatingpoint)));
    CUDAcommon::handleerror(cudaMemset(gpu_Uvec, 0.0, sizeof (floatingpoint)));
#else
    floatingpoint *gpu_Uvec = CUDAcommon::getCUDAvars().gpu_energy;
    /*floatingpoint *gpu_Uvec;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (floatingpoint)));
    CUDAcommon::handleerror(cudaMemsetAsync(gpu_Uvec, 0, sizeof (floatingpoint),streamF));*/
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaMemset(CUDAcommon::cudavars.gpu_energyvec, 0, bntaddvec2.at(0) * sizeof
            (floatingpoint)));
#endif
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
/*    auto gU_tot = CUDAcommon::getCUDAvars().gpu_energy;
    setenergytozero << < 1, 1, 0, streamF >> > (gU_tot);*/
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"ForceFieldManager.cu",
//                            "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    floatingpoint cuda_lambda[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));

#endif
    short count = 0;
    CUDAcommon::tmin.computeenergycalls++;
/*    if(areEqual(d,0.0))
    	CUDAcommon::tmin.computeenerycallszero++;
    else
	    CUDAcommon::tmin.computeenerycallsnonzero++;*/
    for (auto &ff : _forceFields) {
        tbegin = chrono::high_resolution_clock::now();
        auto tempEnergy = ff->computeEnergy(coord, stretched);
        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
//        cout<<ff->getName()<<" "<<tempEnergy<<"pN.nm"<<" ";
        if(CUDAcommon::tmin.individualenergies.size() == _forceFields.size())
            CUDAcommon::tmin.individualenergies[count] += elapsed_energy.count();
        else
            CUDAcommon::tmin.individualenergies.push_back(elapsed_energy.count());

		    if(!stretched){
			    if(CUDAcommon::tmin.individualenergieszero.size() == _forceFields.size())
				    CUDAcommon::tmin.individualenergieszero[count] += elapsed_energy.count();
			    else
				    CUDAcommon::tmin.individualenergieszero.push_back(elapsed_energy.count());
		    }
		    else{
			    if(CUDAcommon::tmin.individualenergiesnonzero.size() == _forceFields.size())
				    CUDAcommon::tmin.individualenergiesnonzero[count] += elapsed_energy
				    		.count();
			    else
				    CUDAcommon::tmin.individualenergiesnonzero.push_back(elapsed_energy
				    .count());
		    }

        count++;
//        cout<<ff->getName()<<" "<<tempEnergy<<endl;
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

        if (verbose) cout << ff->getName() << " energy = " << tempEnergy << endl;
        //if energy is infinity, exit with infinity.
        if (tempEnergy <= -1) {

            //if this is the current energy, exit ungracefully
            if (!stretched) {

                cout << "Energy = " << tempEnergy << endl;

                cout
                        << "Energy of system became infinite. Try adjusting minimization parameters."
                        << endl;
                cout << "The culprit was ... " << ff->getName() << endl;

                _culpritForceField = ff;
                return numeric_limits<floatingpoint>::infinity();
            }
                //if this is a minimization try, just return infinity
            else {
                //cout<<"Returning infintie energy "<<ff->getName()<<endl;
                return numeric_limits<floatingpoint>::infinity();}
        }
        else energy += tempEnergy;
#ifdef SERIAL_CUDACROSSCHECK
        cudaDeviceSynchronize();
        resetfloatingpointvariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
        addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint)
                                                                   , streamF>>>
        (CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                gpu_Uvec);
        floatingpoint cuda_energyvec[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
        std::cout<<ff->getName()<<" Energy. CUDA "<<cuda_energyvec[0]<<" SERL "
                ""<<energy<<endl;
#endif
    }
//    cout<<endl;
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
    //Add energies
#ifdef CUDAACCL
//    std::cout<<"Total nint "<<bntaddvec2.at(0)<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    //Synchronize streams
    for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "computeEnergy",
                                    "ForceFieldManager.cu");
        }
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run2.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run2.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
//    std::cout<<"CUDA energy total nint "<<CUDAcommon::cudavars.offset_E<<endl;
    /*vector<floatingpoint> ones;
    for(int i = 0;i<8192;i++)
        ones.push_back(1.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, ones
                                                        .data() ,
                                                bntaddvec2.at(0) * sizeof
                                                        (floatingpoint),
                                                cudaMemcpyHostToDevice,streamF));*/
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 1,
                                            bntaddvec2.at(0) * sizeof
            (floatingpoint),streamF));
    cudaDeviceSynchronize();*/
    resetfloatingpointvariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
    addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint)
                    , streamF>>>(CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                        gpu_Uvec);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#ifdef DETAILEDOUTPUT_ENERGY
    cudaDeviceSynchronize();
    floatingpoint cuda_energyvec[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"vector energy addition CUDA "<<cuda_energyvec[0]<<" SERL "<<energy<<endl;
#endif
//    CUDAcommon::handleerror(cudaFree(CUDAcommon::cudavars.gpu_energyvec));
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run3(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run3.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run3.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run3.count();
//    std::cout<<"Time total computeEnergy (s) CUDA "<<CUDAcommon::cudatime
//            .TcomputeE<<" SERL "<<CUDAcommon::serltime.TcomputeE<<" factor "
//                     ""<<CUDAcommon::serltime.TcomputeE/CUDAcommon::cudatime.TcomputeE<<endl;
//    std::cout<<"Time split computeEnergy (s) CUDA ";
//    for(auto x:CUDAcommon::cudatime.TveccomputeE)
//        std::cout<<x<<" ";
//    std::cout<<endl;
//    std::cout<<"Time split computeEnergy (s) SERL ";
//    for(auto x:CUDAcommon::serltime.TveccomputeE)
//        std::cout<<x<<" ";
//    std::cout<<endl;
#endif

    return energy;
    
}




EnergyReport ForceFieldManager::computeEnergyHRMD(floatingpoint *coord) const {
    EnergyReport result;
    result.total = 0.0;
    for (auto &ff : _forceFields) {
        auto tempEnergy = ff->computeEnergy(coord);
        // convert to units of kT
        tempEnergy = tempEnergy / kT; // TODO: Energy unit conversion might happen outside
        result.individual.push_back({ff->getName(), tempEnergy});
        result.total += tempEnergy;
    }
    return result;
}




template floatingpoint ForceFieldManager::computeEnergy< false >(floatingpoint *, bool) const;
template floatingpoint ForceFieldManager::computeEnergy< true >(floatingpoint *, bool) const;

void ForceFieldManager::computeForces(floatingpoint *coord, floatingpoint *f) {
    //reset to zero
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    CUDAcommon::cudatime.TcomputeF = 0.0;
    CUDAcommon::cudatime.TveccomputeF.clear();
    CUDAcommon::serltime.TcomputeF = 0.0;
    CUDAcommon::serltime.TveccomputeF.clear();
    tbegin = chrono::high_resolution_clock::now();
#endif
    //@{
    for (int i = 0; i < Bead::getDbData().forces.size_raw(); i++)
        f[i] = 0.0;
    //@}
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::serltime.TcomputeF += elapsed_run.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    CUDAvars cvars = CUDAcommon::getCUDAvars();
    if (cross_checkclass::Aux)
        resetForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, streamF >> >
                                                                      (cvars.gpu_forceAux, gpu_nint);
    else
        resetForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, streamF >> >
                                                                      (cvars.gpu_force, gpu_nint);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));

    CUDAcommon::handleerror(cudaGetLastError(), "resetForcesCUDA", "ForceFieldManager.cu");
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run2.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
    //recompute
//    floatingpoint *F_i = new floatingpoint[CGMethod::N];
    short count = 0;
    CUDAcommon::tmin.computeforcescalls++;
    for (auto &ff : _forceFields) {
        tbegin = chrono::high_resolution_clock::now();
        ff->computeForces(coord, f);
        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
        if(CUDAcommon::tmin.individualforces.size() == _forceFields.size())
            CUDAcommon::tmin.individualforces[count]+= elapsed_energy.count();
        else
            CUDAcommon::tmin.individualforces.push_back(elapsed_energy.count());
        count++;
/*        for(int i=0;i<CGMethod::N;i++){
            if(isnan(f[i])||isinf(f[i])){
                cout<<"Culprit ForceField "<<ff->getName()<<endl;
            }
        }*/
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

//        if(cross_checkclass::Aux)
//            CUDAcommon::handleerror(
//                cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_forceAux, 3 * Bead::getBeads().size() * sizeof
//                                   (floatingpoint),
//                           cudaMemcpyDeviceToHost));
//        else
//            CUDAcommon::handleerror(
//                    cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() * sizeof
//                                       (floatingpoint),
//                               cudaMemcpyDeviceToHost));
//        floatingpoint fmax = 0.0;
//        int id=0;
//        for (auto iter = 0; iter < Bead::getBeads().size(); iter++) {
//            if(abs(F_i[3 *iter])> fmax) {fmax = abs(F_i[3*iter]);id = iter;}
//            if(abs(F_i[3 *iter +1])> fmax) {fmax = abs(F_i[3*iter +1]);id = iter;}
//            if(abs(F_i[3 *iter +2])> fmax) {fmax = abs(F_i[3*iter +2]);id = iter;}
////            std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
//        }
//        std::cout <<"Fmax "<< id<<" "<<fmax<<" "<<F_i[3 * id] << " " << F_i[3 * id + 1] << " " << F_i[3 * id + 2] <<
//                                                                                                                 endl;
    }
//    delete F_i;
}

void ForceFieldManager::computeLoadForces() {

    //reset
    for (auto b: Bead::getBeads()) {
        std::fill(b->loadForcesM.begin(), b->loadForcesM.end(), 0.0);
        std::fill(b->loadForcesP.begin(), b->loadForcesP.end(), 0.0);
//        b->loadForcesP.clear();
//        b->loadForcesM.clear();
    }

    for (auto &f : _forceFields)
        f->computeLoadForces();

    //reset lfi as well
    for (auto b: Bead::getBeads()) {
        b->lfip = 0;
        b->lfim = 0;
    }
}
void ForceFieldManager::computeLoadForce(Cylinder* c, ForceField::LoadForceEnd end) const {
    auto  b          = (end == ForceField::LoadForceEnd::Plus ? c->getSecondBead() : c->getFirstBead());
    auto& loadForces = (end == ForceField::LoadForceEnd::Plus ? b->loadForcesP     : b->loadForcesM   );
    auto& lfi        = (end == ForceField::LoadForceEnd::Plus ? b->lfip            : b->lfim          );

    // reset
    std::fill(loadForces.begin(), loadForces.end(), 0.0);

    for(auto& f : _forceFields) f->computeLoadForce(c, end);

    // reset lfi
    lfi = 0;
}

void ForceFieldManager::printculprit(floatingpoint* force){

    /*cout<<"Printing cylinder data overall"<<endl;
    if(true) {

        cylinder *cylindervec = CUDAcommon::serlvars.cylindervec;
        Cylinder **Cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
        CCylinder **ccylindervec = CUDAcommon::serlvars.ccylindervec;
        floatingpoint *coord = CUDAcommon::serlvars.coord;
        std::cout << "check revectorized cylinders" << endl;
        std::cout << "3 Total Cylinders " << Cylinder::getCylinders().size() << " Beads "
                  << Bead::getBeads().size() << "maxcindex " << Cylinder::getmaxcindex() <<
                  endl;
        for (auto cyl:Cylinder::getCylinders()) {
            int i = cyl->_dcIndex;
            int id1 = cylindervec[i].ID;
            int id2 = Cylinderpointervec[i]->getID();
            int id3 = ccylindervec[i]->getCylinder()->getID();
            if (id1 != id2 || id2 != id3 || id3 != id1)
                std::cout << id1 << " " << id2 << " " << id3 << endl;
            auto b1 = cyl->getFirstBead();
            auto b2 = cyl->getSecondBead();
            long idx1 = b1->_dbIndex;
            long idx2 = b2->_dbIndex;
            cylinder c = cylindervec[i];
            std::cout << "bindices for cyl with ID " << cyl->getID() << " cindex " << i <<
                      " are " << idx1 << " " << idx2 << " " << c.bindices[0] << " "
                      << c.bindices[1] <<" coords ";
            std::cout << coord[3 * idx1] << " " << coord[3 * idx1 + 1] << " "
                      << coord[3 * idx1 + 2] << " " << coord[3 * idx2] << " "
                      << coord[3 * idx2 + 1] << " " << coord[3 * idx2 + 2] <<" forces ";
            std::cout << force[3 * idx1] << " " << force[3 * idx1 + 1] << " "
                      << force[3 * idx1 + 2] << " " << force[3 * idx2] << " "
                      << force[3 * idx2 + 1] << " " << force[3 * idx2 + 2] << endl;
            if (c.bindices[0] != idx1 || c.bindices[1] != idx2) {

                std::cout << "Bead " << b1->coordinate[0] << " " << b1->coordinate[1]
                          << " " << b1->coordinate[2] << " " << " " << b2->coordinate[0]
                          << " " << b2->coordinate[1] << " " << b2->coordinate[2]
                          << " idx " << b1->_dbIndex << " " << b2->_dbIndex << "ID "
                                                                               ""
                          << b1->getID() << " " << b2->getID() << endl;

                std::cout << coord[3 * idx1] << " " << coord[3 * idx1 + 1] << " "
                          << coord[3 * idx1 + 2] << " " << coord[3 * idx2] << " "
                          << coord[3 * idx2 + 1] << " " << coord[3 * idx2 + 2] << endl;
                std::cout << force[3 * idx1] << " " << force[3 * idx1 + 1] << " "
                          << force[3 * idx1 + 2] << " " << force[3 * idx2] << " "
                          << force[3 * idx2 + 1] << " " << force[3 * idx2 + 2] << endl;
                exit(EXIT_FAILURE);
            }
        }

    cout<<"-------DONE------"<<endl;
    }*/

	//get the culprit in output
	_culpritForceField->whoIsCulprit();

	exit(EXIT_FAILURE);
}

#ifdef CUDAACCL

void ForceFieldManager::CUDAcopyForces(cudaStream_t stream, floatingpoint *fprev, floatingpoint *f) {


//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_forceAux));
//    floatingpoint* gpu_forceAux;
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, CGMethod::N * sizeof(floatingpoint)));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.gpu_forceAux=gpu_forceAux;
//    CUDAcommon::cudavars=cvars;

//    std::cout<<"Copyforces Number of Blocks: "<<blocksnthreads[0]<<endl;
//    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;
    copyForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, stream >> >
                                                                 (f, fprev, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "copyForcesCUDA", "ForceFieldManager.cu");
}

#endif

void ForceFieldManager::assignallforcemags() {

    for (auto &ff : _forceFields)
        ff->assignforcemags();
}


void ForceFieldManager::computeHessian(floatingpoint *coord, floatingpoint *f, int total_DOF, float delta){
    // store the minimization time and initialize the matrix
    tauVector.push_back(tau());
    
    chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
    
    vector<vector<floatingpoint> > hessianMatrix(total_DOF, vector<floatingpoint>(total_DOF));
    
    vector<Triplet> tripletList;
    
    
    // loop through all the coordinates
    for(auto i = 0; i < total_DOF; i++){

        // create new vectors for the foorces and coordinates
        vector<floatingpoint> forces_copy_p(total_DOF);
        vector<floatingpoint> coord_copy_p(total_DOF);
        
        vector<floatingpoint> forces_copy_m(total_DOF);
        vector<floatingpoint> coord_copy_m(total_DOF);

        // copy coordinates to new vector
        for(auto l = 0; l< coord_copy_p.size(); l++){
            coord_copy_p[l] = coord[l];
            coord_copy_m[l] = coord[l];
        }

        // perturb the coordinate i
        coord_copy_p[i] += delta;
        coord_copy_m[i] -= delta;

        // calculate the new forces based on perturbation
        computeForces(coord_copy_p.data(), forces_copy_p.data());
        computeForces(coord_copy_m.data(), forces_copy_m.data());

        for(auto j = 0; j < total_DOF; j++){

            // store the derivative of the force on each coordinate j
            float h_i_j = -(forces_copy_p[j] - forces_copy_m[j]) / (2 * delta);
            hessianMatrix[i][j] = h_i_j;
            tripletList.push_back(Triplet(i,j,h_i_j));

        }

    }
    
    // create symmetrized sparse matrix object
    Eigen::SparseMatrix<double> hessMat(total_DOF, total_DOF), hessMatSym;
    hessMat.setFromTriplets(tripletList.begin(), tripletList.end());
    hessMatSym = 0.5*(Eigen::SparseMatrix<double>(hessMat.transpose()) + hessMat);
    
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_vecmat(t1 - t0);
    
    bool denseEstimationBool = SysParams::Mechanics().denseEstimationBool;
    
    if(denseEstimationBool){
    
        Eigen::MatrixXd denseHessMatSym;
        denseHessMatSym = Eigen::MatrixXd(hessMatSym);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> denseESolver(denseHessMatSym);
        evalues = denseESolver.eigenvalues().real();
        evectors = denseESolver.eigenvectors().real();
    
    }else{
        
        Spectra::SparseSymShiftSolve<double> op(hessMatSym);
        //Spectra::SparseSymMatProd<double> op(hessMatSym);
        int numEigs = total_DOF - 1;
        Spectra::SymEigsShiftSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymShiftSolve<double>> eigs(&op, numEigs, total_DOF, 10000);
        //Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double>> eigs(&op, numEigs, numEigs+1);
        /*
        if(evectors.size()!=0){
            //const Eigen::Matrix<double, Eigen::Dynamic, 1> init_vec = evectors.col(0).real();
            const Eigen::Matrix<double, Eigen::Dynamic, 1> init_vec = evectors.real().rowwise().sum();
            if(init_vec.size() == total_DOF){
                const double * arg = init_vec.data();
                eigs.init(arg);
                //eigs.init();
            }else{
                eigs.init();
            };
        }else{
            eigs.init();
        }*/
        
        eigs.init();
        int nconv = eigs.compute();
        evalues = eigs.eigenvalues();
        //columns of evectors matrix are the normalized eigenvectors
        evectors = eigs.eigenvectors(numEigs);
    };
  
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_veceigs(t2 - t1);
    Eigen::VectorXcd IPRI(evectors.cols());
    Eigen::VectorXcd IPRII(evectors.cols());

    // compute participation ratios
    for(auto i = 0; i<evectors.cols(); i++){
        floatingpoint RI = 0.0;
        Eigen::VectorXd col = evectors.col(i).cwiseAbs2();
        RI = pow(col.norm(),2);
        floatingpoint RII = 0.0;
        for(auto j = 0; j < evectors.rows()/3; j++){
            floatingpoint temp = 0.0;
            for(auto k =0; k< 3; k ++){
                temp += pow(evectors(3*j+k,i).real(),2);
            }
            RII += pow(temp,2);
        }
        IPRI(i) = RI;
        IPRII(i) = RII;
    }
    
    
    chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_vecPR(t3 - t2);
    
    cout<<"DOF is "<<total_DOF<<endl;
    std::cout<<"Matrix time "<<elapsed_vecmat.count()<<endl;
    std::cout<<"Compute time "<<elapsed_veceigs.count()<<endl;
    std::cout<<"PR time "<<elapsed_vecPR.count()<<endl;


    // store the full matrix in list
    hessianVector.push_back(hessianMatrix);
    
    // store the eigenvalues in list
    IPRIVector.push_back(IPRI);
    IPRIIVector.push_back(IPRII);
    evaluesVector.push_back(evalues);


}
