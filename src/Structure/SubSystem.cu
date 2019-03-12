
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
#include "MathFunctions.h"
#include "BoundaryElement.h"
#include "BoundaryElementImpl.h"
#include <vector>
#include "dist_driver.h"
#include "dist_coords.h"
#include "dist_common.h"
#include "Cylinder.h"
using namespace mathfunc;
void SubSystem::resetNeighborLists() {
#ifdef CUDAACCL_NL
    coord = new double[CGMethod::N];
                coord_com = new double[3 * Cylinder::getCylinders().size()];
                beadSet = new int[2 * Cylinder::getCylinders().size()];
                cylID = new int[Cylinder::getCylinders().size()];
                filID = new int[Cylinder::getCylinders().size()];
                filType = new int[Cylinder::getCylinders().size()];
                cmpID = new unsigned int[Cylinder::getCylinders().size()];
                fvecpos = new int[Cylinder::getCylinders().size()];
        //        cylstate= new bool[Cylinder::getCylinders().size()];
        //        cylvecpospercmp = new int[2 * _compartmentGrid->getCompartments().size()];

                if(SysParams::Chemistry().numFilaments > 1) {
                    cout << "CUDA NL cannot handle more than one type of filaments." << endl;
                    exit(EXIT_FAILURE);
                }
                int numBindingSites = SysParams::Chemistry().bindingSites[0].size();
                if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                    cmon_state_brancher = new int[ numBindingSites * Cylinder::getCylinders().size()];
                if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                    cmon_state_linker = new int[numBindingSites * Cylinder::getCylinders().size()];
                if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                    cmon_state_motor = new int[numBindingSites * Cylinder::getCylinders().size()];

                int i = 0; //int cID = 0;
                for(auto b: Bead::getBeads()) {
                    //flatten indices
                    int index = 3 * i;
                    coord[index] = b->coordinate[0];
                    coord[index + 1] = b->coordinate[1];
                    coord[index + 2] = b->coordinate[2];
                    i++;
                }
                i = 0; //int countcyl = 0;
        //        for(auto C : _compartmentGrid->getCompartments()) {
        //            int iter = 0;
        //            for(auto nC : C->getNeighbours()) {
        //                cmp_neighbors[nneighbors * cID + iter] = GController::getCompartmentID(nC->coordinates());
        //                        iter++;
        //            }
        //            for(auto k = iter; k < nneighbors; k++ )
        //                cmp_neighbors[nneighbors * cID + k] = -1;

        //            cylvecpospercmp[2 * cID] = countcyl;
        //            countcyl += C->getCylinders().size();
        //            cylvecpospercmp[2 * cID + 1] = countcyl;

        //            if(C->getCylinders().size()>maxnumCyl)
        //                maxnumCyl = C->getCylinders().size();
        //            cID++;
        //            for(auto c:C->getCylinders()){
                 for(auto c:Cylinder::getCylinders()){
                            //flatten indices
                            int index = 3 * i;
                            coord_com[index] = c->coordinate[0];
                            coord_com[index + 1] = c->coordinate[1];
                            coord_com[index + 2] = c->coordinate[2];

                        beadSet[2 * i] = c->getFirstBead()->_dbIndex;
                        beadSet[2 * i + 1] = c->getSecondBead()->_dbIndex;
                        cylID[i] = c->getID();
                        c->_dcIndex = i;
                        fvecpos[i] = c->getPosition();
                        auto fil = dynamic_cast<Filament*>(c->getParent());
                        filID[i] =  fil->getID();
                        cmpID[i] = GController::getCompartmentID(c->getCompartment()->coordinates());
                        filType[i] = fil->getType();
        //                cylstate[i] = c->isFullLength();
                        int j = 0;
                        for(auto it2 = SysParams::Chemistry().bindingSites[fil->getType()].begin();
                            it2 != SysParams::Chemistry().bindingSites[fil->getType()].end(); it2++) {
                            if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                                cmon_state_brancher[numBindingSites * i + j ] = c->getCCylinder()->getCMonomer(*it2)
                                        ->speciesBound(SysParams::Chemistry().brancherBoundIndex[fil->getType()])->getN();
                            if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                                cmon_state_linker[numBindingSites * i + j ] = c->getCCylinder()->getCMonomer(*it2)
                                        ->speciesBound(SysParams::Chemistry().linkerBoundIndex[fil->getType()])->getN();
                            if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                                cmon_state_motor[numBindingSites * i + j ] = c->getCCylinder()->getCMonomer(*it2)
                                        ->speciesBound(SysParams::Chemistry().motorBoundIndex[fil->getType()])->getN();
                            j++;
        //                    for(auto k = 0; k< SysParams::Chemistry().numBoundSpecies[0]; k ++) {
        //                        cmon_state[SysParams::Chemistry().bindingSites[fil->getType()].size() * SysParams::Chemistry
        //                                ().numBoundSpecies[0] * i + j] =
        //                                c->getCCylinder()->getCMonomer(*it2)->speciesBound(k)->getN();
        //                        j++;
        //                    }
                        }
                        i++;
                    }
        //        }//Compartment
                //CUDAMALLOC
                //@{
        //        size_t free, total;
        //        CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
        //        fprintf(stdout,"\t### Available VRAM : %g Mo/ %g Mo(total)\n\n",
        //                free/1e6, total/1e6);
        //
        //        cudaFree(0);
        //
        //        CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
        //        fprintf(stdout,"\t### Available VRAM : %g Mo/ %g Mo(total)\n\n",
        //                free/1e6, total/1e6);
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_coord, CGMethod::N * sizeof(double)),"cuda data "
                                        "transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_coord_com, 3 * Cylinder::getCylinders().size() * sizeof
                                        (double)),"cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, 2 * Cylinder::getCylinders().size() * sizeof
                                        (int)), "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylID, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_fvecpos, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_filID, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_filType, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmpID, Cylinder::getCylinders().size() * sizeof(unsigned
                                        int)), "cuda data transfer", " SubSystem.h");
        //        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylstate, Cylinder::getCylinders().size() * sizeof(int)),
        //                                "cuda data transfer", " SubSystem.h");

                if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmon_state_brancher, numBindingSites *
                        Cylinder::getCylinders().size() * sizeof(int)), "cuda data transfer", " SubSystem.h");
                if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmon_state_linker, numBindingSites *
                                            Cylinder::getCylinders().size() * sizeof(int)), "cuda data transfer", " SubSystem.h");
                if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmon_state_motor, numBindingSites *
                                            Cylinder::getCylinders().size() * sizeof(int)), "cuda data transfer", " SubSystem.h");

        //        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylvecpospercmp,
        //                                2 * _compartmentGrid->getCompartments().size()), "cuda data transfer", " SubSystem.h");
        //        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmp_neighbors, nneighbors *
        //                                 _compartmentGrid->getCompartments().size() * sizeof(int)), "cuda data transfer",
        //                                " SubSystem.h");
                //@}
                //CUDAMEMCPY
                //@{
                CUDAcommon::handleerror(cudaMemcpy(gpu_coord, coord, CGMethod::N *sizeof(double), cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_coord_com, coord_com, 3 * Cylinder::getCylinders().size() *sizeof
                                                   (double), cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, 2 * Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_cylID, cylID, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_fvecpos, fvecpos, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_filID, filID, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_filType, filType, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_cmpID, cmpID, Cylinder::getCylinders().size() *sizeof(unsigned int),
                                                   cudaMemcpyHostToDevice));
        //        CUDAcommon::handleerror(cudaMemcpy(gpu_cylstate, cylstate, Cylinder::getCylinders().size() *sizeof(int),
        //                                           cudaMemcpyHostToDevice));
                if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                CUDAcommon::handleerror(cudaMemcpy(gpu_cmon_state_brancher, cmon_state_brancher, numBindingSites *
                                        Cylinder::getCylinders().size() * sizeof(int), cudaMemcpyHostToDevice));
                if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMemcpy(gpu_cmon_state_linker, cmon_state_linker, numBindingSites *
                                        Cylinder::getCylinders().size() *sizeof(int), cudaMemcpyHostToDevice));
                if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMemcpy(gpu_cmon_state_motor, cmon_state_motor, numBindingSites *
                                        Cylinder::getCylinders().size() *sizeof(int), cudaMemcpyHostToDevice));

        //        CUDAcommon::handleerror(cudaMemcpy(gpu_cylvecpospercmp, cylvecpospercmp,
        //                                           2 * _compartmentGrid->getCompartments().size(), cudaMemcpyHostToDevice));
                CylCylNLvars cylcylnlvars;
                cylcylnlvars.gpu_coord = gpu_coord;
                cylcylnlvars.gpu_coord_com = gpu_coord_com;
                cylcylnlvars.gpu_beadSet = gpu_beadSet;
                cylcylnlvars.gpu_cylID = gpu_cylID;
                cylcylnlvars.gpu_fvecpos = gpu_fvecpos;
                cylcylnlvars.gpu_filID = gpu_filID;
                cylcylnlvars.gpu_filType = gpu_filType;
                cylcylnlvars.gpu_cmpID = gpu_cmpID;
        //        cylcylnlvars.gpu_cylstate = gpu_cylstate;
                cylcylnlvars.gpu_cmon_state_brancher = gpu_cmon_state_brancher;
                cylcylnlvars.gpu_cmon_state_linker = gpu_cmon_state_linker;
                cylcylnlvars.gpu_cmon_state_motor = gpu_cmon_state_motor;
        //        cylcylnlvars.gpu_cylvecpospercmp = gpu_cylvecpospercmp;

                CUDAcommon::cylcylnlvars = cylcylnlvars;
        //        CUDAcommon::handleerror(cudaMemcpy(gpu_cmp_neighbors, cmp_neighbors, nneighbors * _compartmentGrid
        //                                           ->getCompartments().size() *sizeof(int),
        //                                           cudaMemcpyHostToDevice));
                //@}
#endif
    //@{ check begins
    /*cylinder* cylindervec  = CUDAcommon::serlvars.cylindervec;
    Cylinder** Cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
    CCylinder** ccylindervec = CUDAcommon::serlvars.ccylindervec;
    double* coord = CUDAcommon::serlvars.coord;
    for(auto cyl:Cylinder::getCylinders()){
        int i = cyl->_dcIndex;
        int id1 = cylindervec[i].ID;
        int id2 = Cylinderpointervec[i]->getID();
        int id3 = ccylindervec[i]->getCylinder()->getID();
        if(id1 != id2 || id2 != id3 || id3 != id1)
            std::cout<<id1<<" "<<id2<<" "<<id3<<endl;
        auto b1 = cyl->getFirstBead();
        auto b2 = cyl->getSecondBead();
        long idx1 = b1->_dbIndex;
        long idx2 = b2->_dbIndex;
        cylinder c = cylindervec[i];
        std::cout << "4 bindices for cyl with ID "<<cyl->getID()<<" cindex " << i <<
                  " are "<< idx1 << " " << idx2 << " " << c.bindices[0] << " " << c.bindices[1] << endl;
        if(c.bindices[0] != idx1 || c.bindices[1] != idx2) {

            std::cout << "Bead " << b1->coordinate[0] << " " << b1->coordinate[1] << " "
                    "" << b1->coordinate[2] << " " << " " << b2->coordinate[0] << " "
                              "" << b2->coordinate[1] << " " << b2->coordinate[2] << " idx "
                      << b1->_dbIndex << " "
                              "" << b2->_dbIndex << endl;

            std::cout << coord[3 * idx1] << " " << coord[3 * idx1 + 1] << " "
                      << coord[3 * idx1 + 2] << " "
                              "" << coord[3 * idx2] << " " << coord[3 * idx2 + 1] << " "
                      << coord[3 * idx2 + 2] << endl;
        }

    }*/
    //check ends
    chrono::high_resolution_clock::time_point mins, mine;
    mins = chrono::high_resolution_clock::now();

#ifdef HYBRID_NLSTENCILLIST
    _HneighborList->reset();
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_H(mine - mins);
    std::cout<<"H NLSTEN reset time "<<elapsed_H.count()<<endl;
    mins = chrono::high_resolution_clock::now();
    for (auto nlist : __bneighborLists.getElements())
        nlist->reset();
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_B(mine - mins);
    std::cout<<"H NLSTEN B reset time "<<elapsed_B.count()<<endl;

#elif defined(NLORIGINAL) || defined(NLSTENCILLIST)
#ifndef HYBRID_NLSTENCILLIST
    for (auto nl: _neighborLists.getElements())
            nl->reset();
#endif
#endif
    /*mins = chrono::high_resolution_clock::now();
    for (auto nl: _neighborLists.getElements())
        nl->reset();
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_NL(mine - mins);
    std::cout<<"NLSTEN reset time "<<elapsed_NL.count()<<endl;*/

}
void SubSystem::updateBindingManagers() {
#ifdef CUDAACCL_NL
    if(SysParams::Chemistry().numFilaments > 1) {
        cout << "CUDA Binding Manager cannot handle more than one type of filaments." << endl;
        exit(EXIT_FAILURE);
    }
    initializebindingsitesearchCUDA();

    if(CUDAcommon::getCUDAvars().conservestreams)
        numbindmgrs = 0;
    //Calculate binding sites in CUDA
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
        //Linker
        if ((lManager = dynamic_cast<LinkerBindingManager *>(manager.get()))) {
            //calculate all binding Sites.
            getallpossiblelinkerbindingsitesCUDA(lManager, cylcylnlvars.gpu_cmon_state_linker);
        }
        //Motor
        else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
            //calculate all binding Sites.
            getallpossiblemotorbindingsitesCUDA(mManager, cylcylnlvars
                    .gpu_cmon_state_motor);
            }
        //Brancher
        else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {
            //calculate all binding Sites.
            getallpossiblebrancherbindingsitesCUDA(bManager, cylcylnlvars
                    .gpu_cmon_state_brancher);

            }
    }

    //Free vars
    terminatebindingsitesearchCUDA();
    //Assign to respective possible bindings.
    assigntorespectivebindingmanagersCUDA();

//    for(auto gpb:gpu_possibleBindings_vec)
//        CUDAcommon::handleerror(cudaFree(gpb),"cudaFree","SubSystem.cu");
//    for(auto pb:possibleBindings_vec)
//        CUDAcommon::handleerror(cudaFreeHost(pb), "cudaFree", "SubSystem.cu");
//    for(auto np:numpairs_vec)
//        CUDAcommon::handleerror(cudaFreeHost(np),"cudaFree","SubSystem.cu");

    //cudaFree
    endresetCUDA();
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST)
    //vectorize
    SysParams::MParams.speciesboundvec.clear();
    int cidx = 0;
    vector<int> ncylvec(SysParams::CParams.numFilaments);// Number of cylinders
    // corresponding to each filament type.
//    vector<int> bspeciesoffsetvec(SysParams::CParams.numFilaments);
    auto cylvec = Cylinder::getCylinders();
    int ncyl = cylvec.size();
    delete [] cylsqmagnitudevector;
    cylsqmagnitudevector = new double[Cylinder::vectormaxsize];
    unsigned long maxbindingsitespercyl = 0;
    for(auto ftype = 0; ftype < SysParams::CParams.numFilaments; ftype++) {
        maxbindingsitespercyl = max(maxbindingsitespercyl,SysParams::Chemistry()
                .bindingSites[ftype].size());
    }
    long vectorsize = maxbindingsitespercyl * Cylinder::vectormaxsize;
    vector<bool> branchspeciesbound(vectorsize);
    vector<bool> linkerspeciesbound(vectorsize);
    vector<bool> motorspeciesbound(vectorsize);//stores species bound corresponding to each
    // cylinder.

    //set the size of each species bound vector
    fill(branchspeciesbound.begin(),branchspeciesbound.begin()+vectorsize, 0);
    fill(linkerspeciesbound.begin(),linkerspeciesbound.begin()+vectorsize, 0);
    fill(motorspeciesbound.begin(),motorspeciesbound.begin()+vectorsize, 0);

    //fill with appropriate values.
    for (auto cyl: cylvec) {
//        cout<<cyl->_dcIndex<<" "<<cyl->getID()<<endl;
/*        if(cyl->_dcIndex > Cylinder::vectormaxsize)
            std::cout<<"Cindex "<<cyl->_dcIndex<<" greater than vectorsize "
                    ""<<Cylinder::vectormaxsize<<endl;*/
        //cyl->_dcIndex = cidx;
        auto _filamentType = cyl->getType();
        auto x1 = cyl->getFirstBead()->coordinate;
        auto x2 = cyl->getSecondBead()->coordinate;
        vector<double> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
        cylsqmagnitudevector[cyl->_dcIndex] = sqmagnitude(X1X2);
        auto cc = cyl->getCCylinder();
        int idx = 0;
        for (auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
             it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {

            branchspeciesbound[maxbindingsitespercyl * cyl->_dcIndex + idx] =
                    (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN());
            linkerspeciesbound[maxbindingsitespercyl * cyl->_dcIndex + idx] =
                    (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN());
            motorspeciesbound[maxbindingsitespercyl * cyl->_dcIndex + idx] =
                    (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN());
            idx++;
        }
    }
    //@}

    /*for(auto ftype = 0; ftype < SysParams::CParams.numFilaments; ftype++) {
        ncylvec.at(ftype) = cidx;//number of cylinders in each filament.
        bspeciesoffsetvec.at(ftype) = branchspeciesbound.size();
        cidx = 0;
        for (auto cyl: cylvec) {
            //cyl->_dcIndex = cidx;
            auto _filamentType = cyl->getParent()->getType();
            if (_filamentType == ftype) {

                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                vector<double> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
                cylsqmagnitudevector[cyl->_dcIndex] = sqmagnitude(X1X2);
                auto cc = cyl->getCCylinder();
                int idx = 0;
                for (auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
                    branchspeciesbound.push_back (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN());
                    linkerspeciesbound.push_back (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN());
                    motorspeciesbound.push_back (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN());
                    idx++;
                }
                cidx++;
            }
        }
    }
    std::cout<<"max cindex "<<Cylinder::maxcindex<<" removed cylinders "
            ""<<Cylinder::removedcindex.size()<<endl;
    std::cout<<"speciesbound size "<<branchspeciesbound.size()<<endl;*/


    SysParams::MParams.speciesboundvec.push_back(branchspeciesbound);
    SysParams::MParams.speciesboundvec.push_back(linkerspeciesbound);
    SysParams::MParams.speciesboundvec.push_back(motorspeciesbound);
    SysParams::CParams.maxbindingsitespercylinder = maxbindingsitespercyl;
    SysParams::MParams.cylsqmagnitudevector = cylsqmagnitudevector;
//    SysParams::MParams.bsoffsetvec = bspeciesoffsetvec;
    SysParams::MParams.ncylvec = ncylvec;
//    std::cout<<SysParams::Mechanics().speciesboundvec.size()<<endl;
//    std::cout<<motorspeciesbound.size()<<endl;
#endif

    chrono::high_resolution_clock::time_point mins, mine;
    mins = chrono::high_resolution_clock::now();
    //SIMD cylinder update
#ifdef SIMDBINDINGSEARCH2
    minsSIMD = chrono::high_resolution_clock::now();
    for(auto C : _compartmentGrid->getCompartments()) {
        C->SIMDcoordinates();
        C->SIMDcoordinates4linkersearch(1);
        C->SIMDcoordinates4motorsearch(1);
        C->getHybridBindingSearchManager()->resetpossibleBindings();
    }
#endif

    if(!initialize) {
        HybridBindingSearchManager::setdOut();
        initialize = true;
    }

#ifdef SIMDBINDINGSEARCH3
    mineSIMD = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMD2(mineSIMD - minsSIMD);
    SIMDtime += elapsed_runSIMD2.count();
    cout<<"SIMD create time "<<elapsed_runSIMD2.count()<<endl;

    minsSIMD = chrono::high_resolution_clock::now();
    for(auto C : _compartmentGrid->getCompartments()) {
        C->SIMDcoordinates_section();
        C->SIMDcoordinates4linkersearch_section(1);
        C->SIMDcoordinates4motorsearch_section(1);
    }
    mineSIMD = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_SIMDpart(mineSIMD - minsSIMD);
    cout<<"SIMD create time "<<elapsed_SIMDpart.count()<<endl;
#endif

    //PROTOCOL 1 This call calculates Binding pairs according to SIMD protocol V1
    if(false) {
        minsSIMD = chrono::high_resolution_clock::now();
        for (auto C : _compartmentGrid->getCompartments()) {
#ifdef SIMBDINDINGSEARCH

            C->getHybridBindingSearchManager()->updateAllPossibleBindingsstencil();
//        C->getHybridBindingSearchManager()->checkoccupancy(_idvec);
            for (auto &manager : C->getFilamentBindingManagers()) {
#ifdef NLSTENCILLIST
                BranchingManager *bManager;
                if (bManager = dynamic_cast<BranchingManager *>(manager.get()))
                    manager->updateAllPossibleBindingsstencil();
#endif
            }
#else
            for(auto &manager : C->getFilamentBindingManagers()) {
#ifdef NLORIGINAL
                manager->updateAllPossibleBindings();
#endif
#ifdef NLSTENCILLIST
                manager->updateAllPossibleBindingsstencil();
#endif
#if defined(NLORIGINAL) && defined(NLSTENCILLIST)
                manager->crosscheck();
#endif
            }
#endif
        }
        mineSIMD = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runSIMD(mineSIMD - minsSIMD);
        SIMDtime += elapsed_runSIMD.count();
        cout << "SIMD time " << elapsed_runSIMD.count() << endl;
        cout << "find time " << HybridBindingSearchManager::findtime << endl;
    }

    //PRINT
/*    for(auto C : _compartmentGrid->getCompartments()) {
       C->getHybridBindingSearchManager()->printbindingsizes();
    }*/



    //PROTOCOL #2 SIMD V2
/*    for(auto C : _compartmentGrid->getCompartments()) {
        C->getHybridBindingSearchManager()->resetpossibleBindings();
    }*/

    //This call calculates Binding pairs according to SIMD protocol V2
    if(true) {
/*        int totalupn = 0;
        for (auto C : _compartmentGrid->getCompartments()) {
            totalupn += C->getuniquepermuteNeighbours().size();
            cout<<C->getuniquepermuteNeighbours().size()<<" ";
        }
        cout<<endl;
        cout<<"Unique permutation neighbors "<<totalupn<<endl;*/
		#ifdef SIMDBINDINGSEARCH2
        minsSIMD = chrono::high_resolution_clock::now();
        for (auto C : _compartmentGrid->getCompartments()) {

            C->getHybridBindingSearchManager()->updateAllPossibleBindingsstencilSIMDV2();

            /*for(auto &manager : C->getFilamentBindingManagers()) {ad
    #ifdef NLSTENCILLIST
                BranchingManager* bManager;
                if(bManager = dynamic_cast<BranchingManager *>(manager.get()))
                    manager->updateAllPossibleBindingsstencil();
    #endif
            }*/

        }
        //PRINT
/*            for(auto C : _compartmentGrid->getCompartments()) {
                C->getHybridBindingSearchManager()->printbindingsizes();
            }*/
        mineSIMD = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runSIMDV2(mineSIMD - minsSIMD);
        SIMDtimeV2 += elapsed_runSIMDV2.count();
        cout << "SIMDV2 time " << elapsed_runSIMDV2.count() << endl;
        cout << "findV2 time " << HybridBindingSearchManager::findtimeV2 << endl;
        cout << "Append time " << HybridBindingSearchManager::appendtime << endl;
/*        cout << "Time taken to parse SIMD " << HybridBindingSearchManager::SIMDparse1SIMDparse1 << endl;
        cout << "Time taken to merge SIMD " << HybridBindingSearchManager::SIMDparse2 << endl;
        cout << "Time taken to copy to main google map "
                "" << HybridBindingSearchManager::SIMDparse3 << endl;
        cout << "Time taken to update bs "
                "" << HybridBindingSearchManager::SIMDcountbs << endl;*/
        #endif

    }

#ifdef SIMDBINDINGSEARCH3
    for (auto C : _compartmentGrid->getCompartments())
        C->getHybridBindingSearchManager()->resetpossibleBindings();
    minsSIMD = chrono::high_resolution_clock::now();
    HybridBindingSearchManager::findtimeV3 = 0.0;
    HybridBindingSearchManager::SIMDV3appendtime = 0.0;
    for (auto C : _compartmentGrid->getCompartments()) {
        C->getHybridBindingSearchManager()->updateAllPossibleBindingsstencilSIMDV3(0);
    }
    mineSIMD = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMDV3(mineSIMD - minsSIMD);
    cout << "SIMDV3 time " << elapsed_runSIMDV3.count() << endl;
    cout << "findV3 time " << HybridBindingSearchManager::findtimeV3 << endl;
    cout << "Append time " << HybridBindingSearchManager::SIMDV3appendtime << endl;
    cout<<"-------"<<endl;
#endif
    //PROTOCOL #3 This call calculates Binding pairs according to HYBRID protocol
    // (non-SIMD).
#ifdef HYBRID_NLSTENCILLIST
if(false) {
/*    for (auto C : _compartmentGrid->getCompartments()) {
        C->getHybridBindingSearchManager()->resetpossibleBindings();
    }*/

    minsHYBD = chrono::high_resolution_clock::now();
    for (auto C : _compartmentGrid->getCompartments()) {
#ifdef HYBRID_NLSTENCILLIST
        C->getHybridBindingSearchManager()->updateAllPossibleBindingsstencilHYBD();
/*        for (auto &manager : C->getFilamentBindingManagers()) {
#ifdef NLSTENCILLIST
            BranchingManager *bManager;
            if (bManager = dynamic_cast<BranchingManager *>(manager.get()))
                manager->updateAllPossibleBindingsstencil();
#endif
        }*/
    }
#endif
    mineHYBD = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runHYBD(mineHYBD - minsHYBD);
    HYBDtime += elapsed_runHYBD.count();
    cout << "HYBD time " << elapsed_runHYBD.count() << endl;
    cout<<"HYBD map time "<<HybridBindingSearchManager::HYBDappendtime<<endl;
}
#endif

    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_orig(mine - mins);
    std::cout<<"BMgr update time "<<elapsed_orig.count()<<endl;
    //PRINT
/*    for(auto C : _compartmentGrid->getCompartments()) {
        C->getHybridBindingSearchManager()->printbindingsizes();
    }*/
//    exit(EXIT_FAILURE);
}

//OBSOLETE
void SubSystem::vectorizeCylinder() {
    delete [] cylindervec;
    delete [] ccylindervec;
    delete [] cylinderpointervec;
    int Ncyl = Cylinder::getCylinders().size();
    cylindervec = new cylinder[Ncyl];
    ccylindervec = new CCylinder*[Ncyl];
    cylinderpointervec = new Cylinder*[Ncyl];
    //Create cylinder structure
    int i = 0;
    for(auto cyl:Cylinder::getCylinders()){
        //set _dcIndex
        cyl->_dcIndex = i;
        //copy attributes to a structure
        cylindervec[i].filamentID = dynamic_cast<Filament*>(cyl->getParent())->getID();
        cylindervec[i].filamentposition = cyl->getPosition();
        cylindervec[i].bindices[0] = cyl->getFirstBead()->_dbIndex;
        cylindervec[i].bindices[1] = cyl->getSecondBead()->_dbIndex;
        cylindervec[i].cmpID = cyl->getCompartment()->getID();
        cylindervec[i].cindex = i;
        auto coord = cyl->coordinate;
        cylindervec[i].coord[0] = coord[0];
        cylindervec[i].coord[1] = coord[1];
        cylindervec[i].coord[2] = coord[2];
        cylindervec[i].type = cyl->getType();
        cylindervec[i].ID = cyl->getID();
        ccylindervec[i] = cyl->getCCylinder();
        cylinderpointervec[i] = cyl;
        i++;
//        for(int bsc = 0; bsc < nbs; bsc++){
//            double c[3], bead1[3],bead2[3];
//
//            memcpy(bead1, &coord[3*cylindervec[i].bindices[0]], 3 * sizeof(double));
//            memcpy(bead2, &coord[3*cylindervec[i].bindices[1]], 3 * sizeof(double));
//            midPointCoordinate(c,bead1,bead2,bindingsitevec[bsc]);
//            bscoord[12*i+bsc*3] = c[0];
//            bscoord[12*i+bsc*3+1] = c[1];
//            bscoord[12*i+bsc*3+2] = c[2];te<<endl;
//        }
    }
/*    std::cout<<"print for consistency "<<endl;
    for(int idx = 0; idx < Ncyl; idx++) {
        if (cylindervec[idx].cindex != ccylindervec[cylindervec[idx].cindex]->getCylinder()
                ->_dcIndex)
            std::cout << "Fatal mismatch " << cylindervec[idx].cindex << " "
                    ""<<ccylindervec[cylindervec[idx].cindex]->getCylinder()->_dcIndex << endl;
    }*/
    CUDAcommon::serlvars.ccylindervec = ccylindervec;
    CUDAcommon::serlvars.cylindervec = cylindervec;
    CUDAcommon::serlvars.cylinderpointervec = cylinderpointervec;

}

#ifdef CUDAACCL_NL
void SubSystem::initializebindingsitesearchCUDA() {
    //@{ 1. InitializeBSsearch
    //Reset variables
    numpairs_vec.clear();
    possibleBindings_vec.clear();
    gpu_possibleBindings_vec.clear();
//   auto x = CMonomer::_numBSpecies;
//    auto var = SysParams::Chemistry().bmanagerdistances;
    //Malloc params

    //Copy necessary cylinder data to GPU memory
    //@{
    //        if (gpu_params == NULL) {
    int params[3];
    params[0] = SysParams::Chemistry().numBindingSites[0];//filType dependant
    params[1] = 0;//filType dependant
    params[2] = SysParams::Geometry().cylinderNumMon[0];//filType dependant.
    params[3] = Cylinder::getCylinders().size();
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 4 * sizeof(int)),
                            "cuda data transfer", "SubSystem.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params, 4 * sizeof(int),
                                       cudaMemcpyHostToDevice));
//        }
//        if (gpu_bindingSites == NULL) {
    auto bindingSites = SysParams::Chemistry().bindingSites[0];//filType dependant
    int cpu_bindingSites[bindingSites.size()];
    int iii = 0;
    for (auto bs:bindingSites)
    {cpu_bindingSites[iii] = int(bs); iii++;}
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_bindingSites, bindingSites.size() *
                                                                    sizeof(int)), "cuda data transfer", "SubSystem.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_bindingSites, cpu_bindingSites,
                                       bindingSites.size() *  sizeof(int), cudaMemcpyHostToDevice));
//        }
    //@}
}

void SubSystem::getallpossiblelinkerbindingsitesCUDA(LinkerBindingManager* lManager,
                                                     int* cmon_state_linker){
    lManager->assigncudavars();
    cudaStream_t  s;
    if(numbindmgrs + 1 > strvec.size() )
    { cudaStreamCreate(&s); strvec.push_back(s);}
    else
        s = strvec.at(numbindmgrs);
    numbindmgrs++;
//    int *cmon_state_linker = cylcylnlvars.gpu_cmon_state_linker;
    //1. Assign optimal blocks and threads
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
    //2. Calculate binding sites
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
        //call CUDA kernel function
        updateAllPossibleBindingsCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                                  (coord, beadSet, cylID, filID, filType, cmpID, NL, numNLpairs, numpairs,
                                                                                          gpu_params, params2, gpu_possibleBindings, cmon_state_linker,
                                                                                          gpu_bindingSites);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
        //copy binding sites back to CPU
        int *cpu_numpairs, *possibleBindings;
        CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                "Subsystem.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                s),"Copy", "Subsystem.cu");
        CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
        CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
        std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
        numpairs_vec.push_back(cpu_numpairs);
        if(cpu_numpairs[0] > 0) {
            CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
                                                  cudaHostAllocDefault), "Copy", "Subsystem.cu");
            CUDAcommon::handleerror(
                    cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] *
                                                                            sizeof(int), cudaMemcpyDeviceToHost,
                                    s), "Copy", "Subsystem.cu");
            possibleBindings_vec.push_back(possibleBindings);
        }
        gpu_possibleBindings_vec.push_back(gpu_possibleBindings);
//                int cpu_numpairs[1];
//                cudaMemcpy(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost);
//                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
//                int possibleBindings[5 * cpu_numpairs[0]];
//                cudaMemcpy(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
//                           cudaMemcpyDeviceToHost);
    }
    lManager->freecudavars();
    //Free NL numpairs
    CUDAcommon::handleerror(cudaFree(numNLpairs),"cudaFree", "SubSystem.cu");
}

void SubSystem::getallpossiblemotorbindingsitesCUDA(MotorBindingManager* mManager, int*
                                                            cmon_state_motor){
    mManager->assigncudavars();
    cudaStream_t  s;
    if(numbindmgrs + 1 > strvec.size() )
    { cudaStreamCreate(&s); strvec.push_back(s);}
    else
        s = strvec.at(numbindmgrs);
    numbindmgrs++;
//    int *cmon_state_motor = cylcylnlvars.gpu_cmon_state_motor;
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
                                                                                          gpu_params, params2,
                                                                                          gpu_possibleBindings, cmon_state_motor,
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
        numpairs_vec.push_back(cpu_numpairs);
        if(cpu_numpairs[0] > 0) {
            CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
                                                  cudaHostAllocDefault), "Copy", "Subsystem.cu");
            CUDAcommon::handleerror(
                    cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] *
                                                                            sizeof(int), cudaMemcpyDeviceToHost,
                                    s), "Copy", "Subsystem.cu");
            possibleBindings_vec.push_back(possibleBindings);
        }
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
    mManager->freecudavars();
    //Free NL numpairs
    CUDAcommon::handleerror(cudaFree(numNLpairs),"cudaFree", "SubSystem.cu");
}

void SubSystem::getallpossiblebrancherbindingsitesCUDA(BranchingManager* bManager,
                                                       int* cmon_state_brancher) {
    bManager->assigncudavars();
    cudaStream_t  s;
    if(numbindmgrs + 1 > strvec.size() )
    { cudaStreamCreate(&s); strvec.push_back(s);}
    else
        s = strvec.at(numbindmgrs);
    numbindmgrs++;
//    int *cmon_state_brancher = cylcylnlvars.gpu_cmon_state_brancher;
    //2. Assign optimal blocks and threads
    vector<int> blocksnthreads;
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
    int nint = Cylinder::getCylinders().size();
    //int *NL = bManager->getNLCUDA();
    //int *numNLpairs = bManager->getNLsizeCUDA();
    int *numpairs = bManager->getpossiblebindingssizeCUDA();
    double *params2 = bManager->getdistancesCUDA();
    int *zone = bManager->getzoneCUDA();
    //std::cout<<"Total Motor NL size "<<nint<<endl;
    if (nint > 0) {
        //Boundary plane
        auto beList = BoundaryElement::getBoundaryElements();
        int nbe = BoundaryElement::getBoundaryElements().size();
        double *beListplane = new double[4 * nbe];
        double *gpu_beListplane;
        for (int i = 0; i < nbe; i++) {

            if(dynamic_cast<PlaneBoundaryElement*>(beList[i])) {
                double *x = new double[4];
                beList[i]->elementeqn(x);
                beListplane[4 * i] = x[0];
                beListplane[4 * i +1] = x[1];
                beListplane[4 * i +2] = x[2];
                beListplane[4 * i +3] = x[3];
            }
            else{
                cout<<"CUDA cannot handle non-plane type boundaries. Exiting..."<<endl;
                exit(EXIT_FAILURE);
            }
        }
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beListplane, 4 * nbe * sizeof(double)));
        CUDAcommon::handleerror(cudaMemcpy(gpu_beListplane, beListplane, 4 * nbe * sizeof(double), cudaMemcpyHostToDevice));
        //
        int *gpu_possibleBindings;

        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_possibleBindings, SysParams::Chemistry()
                                                                                    .numBindingSites[0] * 3 * nint * sizeof(int)));
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                           updateAllPossibleBindingsCUDA, 0, 0);
        blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreads.push_back(blockSize);
        std::cout << "Brancher blocks and threads " << blocksnthreads.at(0) << " " << blocksnthreads.at(1)
                  << endl;
        resetintvariableCUDA << < 1, 1, 0, s >> > (numpairs);

        updateAllPossibleBindingsBrancherCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                                          (coord, beadSet, cylID, filID, filType, cmpID, numpairs,
                                                                                                  gpu_params, params2, zone, gpu_possibleBindings,gpu_bindingSites,
                                                                                                  cmon_state_brancher, gpu_beListplane);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
        //copy back to CPU
        int *cpu_numpairs, *possibleBindings;
        CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                "Subsystem.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                s),"Copy", "Subsystem.cu");
        CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
        //CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
        std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
        numpairs_vec.push_back(cpu_numpairs);
        if(cpu_numpairs[0] > 0) {
            CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, cpu_numpairs[0] * 3 * sizeof(int),
                                                  cudaHostAllocDefault), "Copy", "Subsystem.cu");
            CUDAcommon::handleerror(
                    cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, cpu_numpairs[0] * 3 *
                                                                            sizeof(int), cudaMemcpyDeviceToHost,
                                    s), "Copy", "Subsystem.cu");
            possibleBindings_vec.push_back(possibleBindings);
        }
        gpu_possibleBindings_vec.push_back(gpu_possibleBindings);
        CUDAcommon::handleerror(cudaFree(gpu_beListplane),"cudaFree", "Subsystem.cu");
    }
    bManager->freecudavars();
}

void SubSystem::terminatebindingsitesearchCUDA(){
    CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "Subsystem.cu");
    CUDAcommon::handleerror(cudaFree(gpu_bindingSites),"cudaFree", "Subsystem.cu");
    //Synchronize streams
    for(auto s:strvec) CUDAcommon::handleerror(cudaStreamSynchronize(s),"stream sync","SubsSystem.cu");
    //Delete sterams
    if(CUDAcommon::getCUDAvars().conservestreams == false) {
        for (auto s:strvec)
            CUDAcommon::handleerror(cudaStreamDestroy(s), "stream destroy", "SubsSystem.cu");
        strvec.clear();
    }
    //clear all possible bindings.
    for(auto c:_compartmentGrid->getCompartments()){
        for(auto &Mgr : c->getFilamentBindingManagers()) {
            Mgr->clearpossibleBindings();
        }
    }
}

void SubSystem::assigntorespectivebindingmanagersCUDA(){
    Compartment* C0 = _compartmentGrid->getCompartments()[0];
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
            if(numpairs[0] > 0)
                count++;
        }
            //MOTORS
        else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[5 * i];
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
            if(numpairs[0] > 0)
                count++;
        }
        else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[3 * i];
                int cIndex = possibleBindings[3 * i +1];
                short cbs = short(possibleBindings[3 * i + 2]);
                auto cylinder = Cylinder::getCylinders()[cIndex];
                Compartment* cmp = GController::getCompartment(cID);
                for(auto &cmanager : cmp->getFilamentBindingManagers()) {
                    if ((bManager = dynamic_cast<BranchingManager *>(cmanager.get()))) {
                        auto t1 = tuple<CCylinder*, short>(cylinder->getCCylinder(), cbs);
                        dynamic_cast<BranchingManager *>(cmanager.get())->appendpossibleBindings(t1);
                    }
                }
            }
//            int n = 0;
//            std::cout<<"-----serial----"<<endl;
//            for(auto C : _compartmentGrid->getCompartments()) {
//                for(auto &bbmanager : C->getFilamentBindingManagers()) {
//                    if ((bManager = dynamic_cast<BranchingManager *>(bbmanager.get()))) {
//                        bbmanager->updateAllPossibleBindings();
//                        n += dynamic_cast<BranchingManager *>(bbmanager.get())->getpossibleBindings().size();
//                    }
//                }
//            }
//            std::cout<<n<<" "<<numpairs_vec[count][0]<<endl;
            if(numpairs[0] > 0)
                count++;
            std::cout<<endl;
        }
    }
}

#endif
CompartmentGrid* SubSystem::_staticgrid;
bool SubSystem::initialize = false;
double SubSystem::SIMDtime  = 0.0;
double SubSystem::SIMDtimeV2  = 0.0;
double SubSystem::HYBDtime  = 0.0;


