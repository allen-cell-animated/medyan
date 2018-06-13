 
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

#ifndef MEDYAN_SubSystem_h
#define MEDYAN_SubSystem_h

#include <vector>
#include <unordered_set>

#include "common.h"

#include "Trackable.h"
#include "Database.h"
#include "Movable.h"
#include "Reactable.h"

#include "NeighborList.h"
#include "DynamicNeighbor.h"

#include "SysParams.h"

#include "CUDAcommon.h"
#include "CGMethod.h"
#include "Filament.h"
#include "Cylinder.h"
#include "CompartmentGrid.h"
#include "Bead.h"
#include "GController.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
//#include "NeighborListImplCUDA.h"

//FORWARD DECLARATIONS
class Boundary;
class Filament;
class Cylinder;
class Linker;
class MotorGhost;
class BranchingPoint;

class CompartmentGrid;

/// Manages all [Movables](@ref Movable) and [Reactables](@ref Reactable). Also holds all
/// [NeighborLists](@ref NeighborList) associated with chemical or mechanical interactions,
/// as well as the CompartmentGrid which contains all chemical structural information, and the
/// system Boundary.

/*! This is a class which handles all changes and information regarding the simulated system.
 *  This class operates as a top manager and provides connections between smaller parts 
 *  of the system. All creation and changes go through this class and will be redirected 
 *  to lower levels. See the databases for more documentation on the explicit creation of
 *  subsystem objects at initialization and during runtime.
 *
 *  This class has functions to add or remove Trackable elements from the system, as well
 *  as update Movable and Reactable instances in the system. It also can update the 
 *  NeighborList container that it holds for the system.
 */
class SubSystem {

public:
    ///Default constructor does nothing
    SubSystem() {}

    ~SubSystem() {}

    /// Add a Trackable to the SubSystem
    template<class T, typename ...Args>
    T *addTrackable(Args &&...args) {

        //create instance
        T *t = new T(forward<Args>(args)...);
        t->addToSubSystem();

        //if movable or reactable, add
        if (t->_movable) addMovable((Movable *) t);

        if (t->_reactable) addReactable((Reactable *) t);

        //if neighbor, add
        if (t->_dneighbor) {
            for (auto nlist : _neighborLists.getElements())
                nlist->addDynamicNeighbor((DynamicNeighbor *) t);
        } else if (t->_neighbor) {
            for (auto nlist : _neighborLists.getElements())
                nlist->addNeighbor((Neighbor *) t);
        }

        return t;
    }

    /// Remove a trackable from the SubSystem
    /// Deleting the actual object should be executed by the actual
    /// callback and/or controlling function.
    template<class T>
    void removeTrackable(T *t) {

        //remove from subsystem
        t->removeFromSubSystem();

        //if movable or reactable, remove
        if (t->_movable) removeMovable((Movable *) t);

        if (t->_reactable) removeReactable((Reactable *) t);

        //if neighbor, remove
        if (t->_dneighbor) {
            for (auto nlist : _neighborLists.getElements())
                nlist->removeDynamicNeighbor((DynamicNeighbor *) t);
        } else if (t->_neighbor) {
            for (auto nlist : _neighborLists.getElements())
                nlist->removeNeighbor((Neighbor *) t);
        }
    }

    //@{
    /// Setter functions for Movable
    void addMovable(Movable *mov) { _movables.insert(mov); }

    void removeMovable(Movable *mov) {
        auto it = _movables.find(mov);
        if (it != _movables.end()) _movables.erase(it);
    }

    //@}
    /// Get all Movable
    const unordered_set<Movable *> &getMovables() { return _movables; }

    //@{
    /// Setter function for Reactable
    void addReactable(Reactable *r) { _reactables.insert(r); }

    void removeReactable(Reactable *r) {
        auto it = _reactables.find(r);
        if (it != _reactables.end()) _reactables.erase(it);
    }

    //@}
    /// Get all Reactable
    const unordered_set<Reactable *> &getReactables() { return _reactables; }

    /// Get the subsystem boundary
    Boundary *getBoundary() { return _boundary; }

    /// Add a boundary to this subsystem
    void addBoundary(Boundary *boundary) { _boundary = boundary; }

    /// Add a neighbor list to the subsystem
    void addNeighborList(NeighborList *nl) { _neighborLists.addElement(nl); }

    /// Reset all neighbor lists in subsystem
    void resetNeighborLists() {
#ifdef CUDAACCL_NL
        //        nvtxRangePushA("NL_Prep_SubSystem");
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
        //        nvtxRangePop();
#endif
//        std::cout<<_neighborLists.getElements().size()<<endl;
        for (auto nl: _neighborLists.getElements())
            nl->reset();
    }

#ifdef CUDAACCL_NL
    void endresetCUDA(){
        for(auto gpb:gpu_possibleBindings_vec)
            CUDAcommon::handleerror(cudaFree(gpb),"cudaFree","SubSystem.cu");
        for(auto pb:possibleBindings_vec)
            CUDAcommon::handleerror(cudaFreeHost(pb), "cudaFree", "SubSystem.cu");
        for(auto np:numpairs_vec)
            CUDAcommon::handleerror(cudaFreeHost(np),"cudaFree","SubSystem.cu");
        auto cylcylnvars = CUDAcommon::getCylCylNLvars();
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_coord),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_coord_com),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_beadSet),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_cylID),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_fvecpos),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_filID),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_filType),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_cmpID),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_cmon_state_linker),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_cmon_state_brancher),"cudaFree","SubSystem.h");
        CUDAcommon::handleerror(cudaFree(cylcylnvars.gpu_cmon_state_motor),"cudaFree","SubSystem.h");
//        Compartment* C0 = _compartmentGrid->getCompartments()[0];

    }
#endif
    //@{
    ///Subsystem energy management
    double getSubSystemEnergy() {return _energy;}
    void setSubSystemEnergy(double energy) {_energy = energy;}
    //@}
    
    //@{
    /// CompartmentGrid management
    void setCompartmentGrid(CompartmentGrid* grid) {_compartmentGrid = grid;}
    CompartmentGrid* getCompartmentGrid() {return _compartmentGrid;}
    //@]
    
    /// Update the binding managers of the system
    void updateBindingManagers();
    
private:

    double _energy = 0; ///< Energy of this subsystem
    Boundary* _boundary; ///< Boundary pointer
    
    unordered_set<Movable*> _movables; ///< All movables in the subsystem
    unordered_set<Reactable*> _reactables; ///< All reactables in the subsystem
        
    Database<NeighborList*> _neighborLists; ///< All neighborlists in the system
        
    CompartmentGrid* _compartmentGrid; ///< The compartment grid
#ifdef CUDAACCL_NL
    double* gpu_coord;
    double* gpu_coord_com;
    int * gpu_beadSet;
    int *gpu_cylID;
    int *gpu_filID;
    int *gpu_cmpID;
//    int *gpu_cylstate;
    int *gpu_cmon_state_brancher;
    int *gpu_cmon_state_linker;
    int *gpu_cmon_state_motor;
//    int *gpu_cylvecpospercmp;
    int *gpu_fvecpos;
    int *gpu_filType;
    double *coord;
    double *coord_com;
    int *beadSet;
    int *cylID;
    int *filID;
    int *filType;
    int *gpu_bindingSites;
    unsigned int *cmpID;
//    bool *cylstate;
    int *cmon_state_brancher;
    int *cmon_state_linker;
    int *cmon_state_motor;
//    int *cylvecpospercmp;//each compartment gets a start and end integer representing the position of first and last
    // cylinder in cmpID;
    int *fvecpos; //position of cylinder in filament vector.
//    int maxnumCyl = 0;
    int *gpu_params = NULL;//used for possiblebindings update
    vector<cudaStream_t > strvec;
    int numbindmgrs = 0;
    vector<int*> numpairs_vec;//Stores number of BindingSite pairs for each binding manager.
    vector<int*> possibleBindings_vec;
    vector<int*> gpu_possibleBindings_vec;
    void initializebindingsitesearchCUDA();
    void getallpossiblelinkerbindingsitesCUDA(LinkerBindingManager* lManager, int* cmon_state_linker);
    void getallpossiblemotorbindingsitesCUDA(MotorBindingManager* mManager, int*
    cmon_state_motor);
    void getallpossiblebrancherbindingsitesCUDA(BranchingManager* bManager, int*
    cmon_state_brancher);
    void terminatebindingsitesearchCUDA();
    void assigntorespectivebindingmanagersCUDA();
#endif
};

#endif
