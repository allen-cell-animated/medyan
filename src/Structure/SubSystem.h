 
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

#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
#include "CGMethod.h"
#include "Filament.h"
#include "Cylinder.h"
#include "CompartmentGrid.h"
#include "Bead.h"
#include "GController.h"
#include "HybridNeighborList.h"
#include "HybridNeighborListImpl.h"

#include <initializer_list>

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
    ///constructor creates and stores an instance of HybridCylinderCylinderNL

    SubSystem() {
#ifdef HYBRID_NLSTENCILLIST
        _HneighborList = new HybridCylinderCylinderNL();
#endif
}

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
#ifdef HYBRID_NLSTENCILLIST
            _HneighborList->addDynamicNeighbor((DynamicNeighbor *) t);
            //Remove boundary and bubble neighbors
            for (auto nlist : __bneighborLists)
                nlist->addDynamicNeighbor((DynamicNeighbor *) t);
#else
            for (auto nlist : _neighborLists)
                nlist->addDynamicNeighbor((DynamicNeighbor *) t);
#endif

        } else if (t->_neighbor) {
            for (auto nlist : _neighborLists)
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
#ifdef HYBRID_NLSTENCILLIST
            _HneighborList->removeDynamicNeighbor((DynamicNeighbor *) t);
            //Remove boundary neighbors
            for (auto nlist : __bneighborLists)
                nlist->removeDynamicNeighbor((DynamicNeighbor *) t);
#else
            for (auto nlist : _neighborLists)
                nlist->removeDynamicNeighbor((DynamicNeighbor *) t);
#endif

        } else if (t->_neighbor) {
            for (auto nlist : _neighborLists)
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
    void addNeighborList(NeighborList *nl) { _neighborLists.push_back(nl); }

    void addBNeighborList(NeighborList *nl) { __bneighborLists.push_back(nl); }

    /// Reset all neighbor lists in subsystem
    void resetNeighborLists();

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
    void setCompartmentGrid(CompartmentGrid* grid) {_compartmentGrid = grid; _staticgrid = _compartmentGrid;}
    CompartmentGrid* getCompartmentGrid() {return _compartmentGrid;}
    //@]
    
    /// Update the binding managers of the system
    void updateBindingManagers();

#ifdef HYBRID_NLSTENCILLIST
    //getter for HneighborList
    HybridCylinderCylinderNL* getHNeighborList(){return _HneighborList;}

    void initializeHNeighborList(){_HneighborList->initializeHybridNeighborList();}
#endif
    static CompartmentGrid* getstaticgrid(){
        return _staticgrid;
    }

    static double HYBDtime;
private:
    double _energy = 0; ///< Energy of this subsystem
    Boundary* _boundary; ///< Boundary pointer
    
    unordered_set<Movable*> _movables; ///< All movables in the subsystem
    unordered_set<Reactable*> _reactables; ///< All reactables in the subsystem
        
    std::vector<NeighborList*> _neighborLists; ///< All neighborlists in the system
    std::vector<NeighborList*> __bneighborLists; ///< Boundary neighborlists in the system.
    // Used only in Hybrid binding Manager cases
#ifdef HYBRID_NLSTENCILLIST
    HybridCylinderCylinderNL* _HneighborList;
#endif
        
    CompartmentGrid* _compartmentGrid; ///< The compartment grid

    //Cylinder vector
    static CompartmentGrid* _staticgrid;
    double* cylsqmagnitudevector = NULL;
    static bool initialize;
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
