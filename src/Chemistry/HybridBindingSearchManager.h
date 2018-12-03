
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
#ifndef MEDYAN_HybridBindingSearchManager_h
#define MEDYAN_HybridBindingSearchManager_h
#ifdef HYBRID_NLSTENCILLIST
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>

#include "common.h"

#include "NeighborListImpl.h"
#include "HybridNeighborListImpl.h"
#include "ReactionBase.h"

#include "SysParams.h"
#include "Rand.h"
#include "BindingManager.h"

//FORWARD DECLARATIONS
class SubSystem;
class ReactionBase;
class CCylinder;
class Compartment;
class Cylinder;
class FilamentBindingManager;

class HybridBindingSearchManager {

    friend class ChemManager;

private:
    Compartment* _compartment;
//    float _rMin; ///< Minimum reaction range
//    float _rMax; ///< Maximum reaction range
//    float _rMinsq; ///< Minimum reaction range squared
//    float _rMaxsq; ///< Maximum reaction range squared
    vector<vector<double>> bindingsites1;
    vector<vector<double>> bindingsites2;
    vector<vector<float>> _rMaxsqvec; //squared maxdistance cutoff
    vector<vector<float>> _rMinsqvec;//squared mindistance cutoff
    vector<vector<short>> _filamentIDvec;//filament ID pairs considered
    static vector<short> HNLIDvec; //Hybrid NL ID to track the total number of
    // neighborilists in the system
    vector<vector<short>> bstateposvec;
    vector<vector<float>> minparamcyl2;
    vector<vector<float>> maxparamcyl2;
    int totaluniquefIDpairs = 0;
    vector<vector<FilamentBindingManager*>> fManagervec;
    //possible bindings at current state. updated according to neighbor list

    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
            _possibleBindings;
    //possible bindings at current state. updated according to neighbor list stencil
    vector<vector<unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>>>
            _possibleBindingsstencilvec;
    vector<vector<unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*,
    short>>>>>_reversepossibleBindingsstencilvec;
    //static neighbor list

    static HybridCylinderCylinderNL* _HneighborList;
public:
    //constructors
     HybridBindingSearchManager(Compartment* compartment);
    ~HybridBindingSearchManager() {}

    //@{
    //@{
    /// Getters for distances
//    float getRMin() {return _rMin;}
//    float getRMax() {return _rMax;}
    //@}

    bool isConsistent();


    void addPossibleBindingsstencil(short idvec[2], CCylinder* cc, short bindingSite);
    void removePossibleBindingsstencil(short idvec[2], CCylinder* cc, short bindingSite);
    ///update all possible binding reactions that could occur using stencil NL
    void updateAllPossibleBindingsstencil();
    vector<tuple<CCylinder*, short>> chooseBindingSitesstencil(short idvec[2]){
        short idx = idvec[0];
        short idx2 = idvec[1];
        int pbsSize = _possibleBindingsstencilvec[idx][idx2].size() ;
        assert((pbsSize!= 0)
               && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");
        int randomIndex = Rand::randInteger(0, pbsSize - 1);
        auto it = _possibleBindingsstencilvec[idx][idx2].begin();

        advance(it, randomIndex);

        return vector<tuple<CCylinder*, short>>{it->first, it->second};

    };
    void setbindingsearchparameter( FilamentBindingManager* fmanager, short bstatepos,
                                    short ftype1 = 0, short ftype2 = 0, float rMax = 0.0, float
                                    rMin = 0.0);

    void addtoHNeighborList();

    vector<Cylinder*> getHNeighbors(Cylinder* c, short HNLID){
        return _HneighborList->getNeighborsstencil(HNLID, c);
    }

    void checkoccupancy(short idvec[2]);

private:


};

#endif
#endif
