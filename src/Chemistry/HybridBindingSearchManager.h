
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

/*#include <deque>
#include <tuple>
#include <thrust/binary_search.h>*/


//using __gnu_cxx::hash;  // or __gnu_cxx::hash, or maybe tr1::hash, depending on your OS
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
    /*struct gdmapstruct{
        gdmap gdmap1;
        gdmap gdmap2;
        gdmap gdmap3;
        gdmap gdmap4;
        gdmap gdmap5;
        gdmap gdmap6;
        gdmap gdmap7;
        gdmap gdmap8;
        gdmap gdmap9;
        gdmap gdmap10;
        gdmap gdmap11;
        gdmap gdmap12;
        gdmap gdmap13;
        gdmap gdmap14;
        gdmap gdmap15;
        gdmap gdmap16;
        gdmap& getelement(int i){
            switch(i){
                case 1: return gdmap1;
                case 2: return gdmap2;
                case 3: return gdmap3;
                case 4: return gdmap4;
                case 5: return gdmap5;
                case 6: return gdmap6;
                case 7: return gdmap7;
                case 8: return gdmap8;
                case 9: return gdmap9;
                case 10: return gdmap10;
                case 11: return gdmap11;
                case 12: return gdmap12;
                case 13: return gdmap13;
                case 14: return gdmap14;
                case 15: return gdmap15;
                case 16: return gdmap16;
            }
        }
    };
    gdmapstruct _gdmapstruct;

    template<typename T>
    class SortOrder
    {
    public:
        SortOrder(const std::vector<T> *_sortArray) : sortArray(_sortArray) {;}

        bool operator()(int lhs, int rhs) const
        {
//            cout<<"pos "<<lhs<<" "<<rhs<<endl;
//            cout<<"arr "<<sortArray[lhs]<<" "<<sortArray[rhs]<<endl;
            return sortArray[lhs] < sortArray[rhs];
        }

    private:
        const std::vector<T> *sortArray;
    };

    template <typename Container>
    struct compare_indirect_index
    {
        const Container& container;
        compare_indirect_index( const Container& container ): container( container ) { }
        bool operator () ( size_t lindex, size_t rindex ) const
        {
            return container[ lindex ] < container[ rindex ];
        }
    };

    bool compareX(int a, int b, uint32_t* data)
    {
        return data[a]<data[b];
    }*/

//    float _rMin; ///< Minimum reaction range
//    float _rMax; ///< Maximum reaction range
//    float _rMinsq; ///< Minimum reaction range squared
//    float _rMaxsq; ///< Maximum reaction range squared

    chrono::high_resolution_clock::time_point minsHYBD, mineHYBD,
            minsfind, minefind, minsmap, minemap;
    Compartment* _compartment;

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
    vector<vector<int>> Nbindingpairs;
    int totaluniquefIDpairs = 0;

    vector<vector<FilamentBindingManager*>> fManagervec;

    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
            _possibleBindings;

/*    possible bindings at current state. updated according to neighbor list stencil
    vector<vector<unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>>>
            _possibleBindingsstencilvec;*/

    vector<vector<unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*,
    short>>>>>_reversepossibleBindingsstencilvec;
    //static neighbor list
    static HybridCylinderCylinderNL* _HneighborList;
    

public:
    //possible bindings at current state. updated according to neighbor list stencil
    vector<vector<unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>>>
            _possibleBindingsstencilvec;


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
    void updateAllPossibleBindingsstencilHYBD();

    vector<tuple<CCylinder*, short>> chooseBindingSitesstencil(short idvec[2]);

    void setbindingsearchparameter(FilamentBindingManager* fmanager, short bstatepos,
                                   short ftype1 = 0, short ftype2 = 0, float rMax = 0.0,
                                   float rMin = 0.0);

    void addtoHNeighborList();

    vector<Cylinder*> getHNeighbors(Cylinder* c, short HNLID){
        return _HneighborList->getNeighborsstencil(HNLID, c);
    }

    void checkoccupancy(short idvec[2]);
    void printbindingsizes(){
        int idx, idx2;
        for(idx = 0; idx<totaluniquefIDpairs; idx++){
            int countbounds = _rMaxsqvec[idx].size();
            for (idx2 = 0; idx2 < countbounds; idx2++) {
                std::cout<<"Hybrid "<<_possibleBindingsstencilvec[idx][idx2].size()<< " SIMD "
                                          ""<<Nbindingpairs[idx][idx2]<<endl;
            }
        }
    }

    void resetpossibleBindings(){

        int idx, idx2;
        for(idx = 0; idx<totaluniquefIDpairs; idx++){
            int countbounds = _rMaxsqvec[idx].size();
            for (idx2 = 0; idx2 < countbounds; idx2++) {

                _possibleBindingsstencilvec[idx][idx2].clear();
                _reversepossibleBindingsstencilvec[idx][idx2].clear();

            }
        }
    }


    static double HYBDtime;
    static double HYBDappendtime;
};
#endif
#endif
