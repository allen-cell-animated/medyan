
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
#ifdef SIMDBINDINGSEARCH
#include <sparsehash/dense_hash_map>
#include "dist_driver.h"

#endif
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
#ifdef SIMDBINDINGSEARCH
using google::dense_hash_map;      // namespace where class lives by default
typedef dense_hash_map<uint32_t , vector<uint32_t>, hash<uint32_t>> gdmap;
typedef dense_hash_map<CCylinder* , vector<uint32_t>, hash<uint32_t>> gdmap2;
typedef tuple<gdmap, gdmap, gdmap, gdmap, gdmap, gdmap, gdmap, gdmap, gdmap, gdmap,
        gdmap, gdmap, gdmap, gdmap, gdmap, gdmap> TUPLEGDMAP;
#endif

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

    chrono::high_resolution_clock::time_point minsSIMD, mineSIMD, minsHYBD, mineHYBD,
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
#ifdef SIMDBINDINGSEARCH


    vector<uint32_t> linker1, linker2;
    vector<uint32_t> motor1, motor2;
    uint Npairs = 0;

    vector<vector<gdmap>> googlepossible;
    vector<vector<gdmap>>  googlereversepossible;
    gdmap googlepossiblel;
    gdmap googlepossiblem;
    gdmap googlereversepossiblel;
    gdmap googlereversepossiblem;
    gdmap googlepossibleORIG;
    gdmap googlereversepossibleORIG;
#endif
    //static neighbor list
    static HybridCylinderCylinderNL* _HneighborList;
    
#ifdef SIMDBINDINGSEARCH
    //SIMD variables
    unsigned mask = (1 << 4) - 1;
    static const dist::tag_simd<dist::simd_avx_par,  float>  t_avx_par;
    static const dist::tag_simd<dist::simd_avx,  float>   t_avx;
    static const dist::tag_simd<dist::simd_no,   float>   t_serial;
    bool initialized = false;

    //D = 1
    static dist::dOut<1U,true> bspairs1self;
    static dist::dOut<1U,false> bspairs1;
    static dist::dOut<1U,false> bspairslinker;
    static dist::dOut<1U,true> bspairslinkerself;
    static dist::dOut<1U,false> bspairsmotor;
    static dist::dOut<1U,true> bspairsmotorself;
    vector<uint32_t> pairslinker;
    vector<uint32_t> pairsmotor;
    vector<bool> filID_fpos_pairL;
    vector<bool> filID_fpos_pairM;
    vector<vector<bool>> pairvaluespecieslinker;
    vector<vector<bool>> pairvaluespeciesmotor;

    template <uint D, bool SELF, bool LinkerorMotor>
    void calculatebspairsLMself(dist::dOut<D, SELF>& bspairs, short idvec[2]);

    template <uint D, bool SELF, bool LinkerorMotor>
    void calculatebspairsLMenclosed(dist::dOut<D, SELF>& bspairs, short idvec[2]);

    template<uint D, bool SELF, bool LinkerorMotor>
    void parseSIMDout(dist::dOut<D,SELF>& bspairsoutS, short idvec[2],
            Compartment* ncmp = NULL);

    template<uint D, bool SELF=true, bool LinkerorMotor>
    void threadwritepairsself(uint startID, uint endID, dist::dOut<D,SELF>& bspairsoutS,
            short idvec[2], gdmap tempmap[], short threadID);

    template<uint D, bool SELF=false, bool LinkerorMotor>
    void threadwritepairs(Compartment* ncmp, dist::dOut<D,SELF>& bspairsoutS,
            gdmap tempmap[], gdmap tempNmap[], uint startID, uint endID,
            short idvec[2], short threadID );

    //V2
    template<uint D, bool SELF=true, bool LinkerorMotor>
    void threadwritepairsselfV2(uint startID, uint endID, dist::dOut<D,SELF>& bspairsoutS,
                              short idvec[2],
                              dense_hash_map<uint32_t , uint, hash<uint32_t>> tempmap[],
                              vector<vector<uint32_t>> valuematrixvec[], short threadID);

    //V3
    template<uint D, bool SELF=true, bool LinkerorMotor>
    void threadwritepairsselfV3(uint startID, uint endID,
            dist::dOut<D, SELF> &bspairsoutS, short idvec[2], short threadID);


    template<uint D, bool SELF=false, bool LinkerorMotor>
    void threadwritepairsV3(Compartment* ncmp,
    dist::dOut<D, SELF> &bspairsoutS, uint startID, uint endID, short idvec[2],
                        short threadID);


    template<uint D, bool LinkerorMotor>
    void threadedsingleparse(uint prev, uint next, short i);

/*    template<uint D, bool LinkerorMotor>
    void threadedsingleparseV2(gdmap tempmap[],
                               vector<uint32_t>& cmpID, vector<uint32_t>& ncmpID,
                               uint prev, uint next, short i);*/

    template<bool LinkerorMotor>
    void mergemaps(short threadID, short offset);

    template<bool LinkerorMotor>
    void findIDsincylinderIDvector(vector<vector<uint>>& outputvector,
            vector<vector<uint32_t>>& cmpIDbs);


    template<uint D, bool LinkerorMotor>
    void singlepassparseSIMDout(short idvec[2]);

    static const short nthreads = 1;
    gdmap vecmapL[2 * nthreads];
    gdmap vecmapM[2 * nthreads];

    deque<gdmap> dqvecmap;

    TUPLEGDMAP tuplegdmap;

    dense_hash_map<uint32_t , uint, hash<uint32_t>> vecmapV2[2 * nthreads];
    dense_hash_map<uint32_t , uint, hash<uint32_t>> vecNmapV2[2 * nthreads];

    vector<vector<uint32_t>> valuematrixvec[2*nthreads];

    template<bool LinkerorMotor>
    vector<uint>& getLinkerorMotor(const short oneortwo){
        if(LinkerorMotor){
            if(oneortwo == 1)
                return linker1;
            else
                return linker2;
        }
        else{
            if(oneortwo == 1)
                return motor1;
            else
                return motor2;
        }
    }

    //getters for vectors.
    template<bool LinkerorMotor>
    vector<bool>& getfilID_fpospairs(){
        if(LinkerorMotor){
            return filID_fpos_pairL;
        }
        else{
            return filID_fpos_pairM;
        }
    }

    template<bool LinkerorMotor>
    vector<uint32_t>& getpairsLinkerorMotor(){
        if(LinkerorMotor){
            return pairslinker;
        }
        else{
            return pairsmotor;
        }
    }

    template<bool LinkerorMotor>
    vector<vector<bool>>& getboundspeciesLinkerorMotor(){
        if(LinkerorMotor){
            return pairvaluespecieslinker;
        }
        else{
            return pairvaluespecieslinker;
        }
    }

    //Getters for vectors of google dense hash map
    template<bool LinkerorMotor>
    gdmap* getvecmapLinkerorMotor(){
        if(LinkerorMotor){
            return vecmapL;
        }
        else{
            return vecmapM;
        }
    }

    //Getters for map
    template<bool LinkerorMotor>
    gdmap& getgooglemapLinkerorMotor(){
        if(LinkerorMotor){
            return googlepossiblel;
        }
        else{
            return googlepossiblem;
        }
    }

    //Getters for reverse map
    template<bool LinkerorMotor>
    gdmap& getgoogleRmapLinkerorMotor(){
        if(LinkerorMotor){
            return googlereversepossiblel;
        }
        else{
            return googlereversepossiblem;
        }
    }

    template<bool LinkerorMotor>
    void copytogooglemap(){
        short partnerID = 0;
        if (getvecmapLinkerorMotor<LinkerorMotor>()[partnerID].size() ) {
            for (const auto &iter : getvecmapLinkerorMotor<LinkerorMotor>()[partnerID]) {
                getgooglemapLinkerorMotor<LinkerorMotor>()[iter.first] = iter.second;
            }
        }
        //reverse map
        if (getvecmapLinkerorMotor<LinkerorMotor>()[partnerID + 1].size()) {
            for (const auto &iter : getvecmapLinkerorMotor<LinkerorMotor>()[partnerID + 1]) {
                getgoogleRmapLinkerorMotor<LinkerorMotor>()[iter.first] = iter.second;
            }
        }
    }

    template <uint D, bool SELF, bool LinkerorMotor>
    void gatherCylinderIDfromcIndex(dist::dOut<D,SELF>& bspairsoutS, int start, int end,
            uint prev_size, Compartment* nCmp = NULL);

    //D = 2
    static dist::dOut<2U,true> bspairs2self;
    static dist::dOut<2U,false> bspairs2;

/*    static dist::dOut<1U,false> bspairs2_D1;
    static dist::dOut<1U,false> bspairs2_D2;
    //D = 3
    dist::dOut<3U,true> bspairs3self;
    dist::dOut<3U,false> bspairs3;*/
    //D = 4
    /*dist::dOut<4U,true> bspairs4self;
    dist::dOut<4U,false> bspairs4;*/

    template <uint D, bool SELF>
    void calculatebspairsself(dist::dOut<D, SELF>& bspairs);
    template <uint D, bool SELF>
    void calculatebspairsenclosed(dist::dOut<D, SELF>& bspairs);

    template <bool LinkerorMotor>
    void calculatebspairsforacylinder(short idvec[2], uint cIndex, CCylinder* cCylinder,
                                      short bindingSite, dist::dOut<1, true>& bspairsself,
                                      dist::dOut<1, false>& bspairs, dist::Coords& coord);

    template <uint D, bool SELF>
    void parsebspairsforacylinder(short idvec[2], dist::dOut<D, SELF>& bspairs,
                                  Compartment* ncmp, uint cIndex1, short bindingsite,
                                  short complimentaryfID);

    static short Totallinkermotor;

    /*template <uint D, bool SELF>
    dist::dOut<D,SELF> bspairs(){
        if(D == 1){
            if(SELF == true) return bspairs1self;
            else return bspairs1;
        }
        else if(D == 2){
            if(SELF == true) return bspairs2self;
            else return bspairs2;
        }
        *//*else if(D == 3){
            if(SELF == true) return bspairs3self;
            else return bspairs3;
        }
        else if(D == 4){
            if(SELF == true) return bspairs4self;
            else return bspairs4;
        }*//*
        else{
            cout<<"Number of distance pairs mentioned in chemistryinput.txt exceeds "
                  "maximum pairs allowed. Exiting.."<<endl;
            exit(EXIT_FAILURE);
        }
    }*/

   void countKeys(vector<uint16_t>& keyvec, vector<unsigned int long>& bounds,
                  vector<int>&  countvec, vector<uint32_t>& bindingskey,
                  uint prev, uint next);

   void countNpairsfound(short idvec[2]);
#endif
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
    void updateAllPossibleBindingsstencilSIMDV2();
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
#ifdef SIMDBINDINGSEARCH
    void checkoccupancySIMD(short idvec[2]);

    void checkoccupancySIMD(short idvec[2], gdmap& map);
#endif
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
#ifdef SIMDBINDINGSEARCH
                googlepossible[idx][idx2].clear();
                googlereversepossible[idx][idx2].clear();
                getfilID_fpospairs<true>().clear();
                getfilID_fpospairs<false>().clear();
                getpairsLinkerorMotor<true>().clear();
                getpairsLinkerorMotor<false>().clear();
#endif
            }
        }
    }
#ifdef SIMDBINDINGSEARCH
    static void setdOut(){
        Totallinkermotor = 2;
        bspairs2self.init_dout(10000, {900.0f, 1600.0f, 30625.0f, 50625.0f});
        bspairs2.init_dout(10000, {900.0f, 1600.0f, 30625.0f, 50625.0f});

/*        bspairs2_D1.init_dout(10000,{900.0f,1600.0f});
        bspairs2_D2.init_dout(10000,{30625.0f, 50625.0f});*/

        // V2
        bspairslinkerself.init_dout(10000,{900.0f,1600.0f});
        bspairslinker.init_dout(10000,{900.0f,1600.0f});
        bspairsmotorself.init_dout(10000,{30625.0f, 50625.0f});
        bspairsmotor.init_dout(10000,{30625.0f, 50625.0f});
    }
#endif

    static double SIMDtime;
    static double HYBDtime;
    static double findtime;
    static double findtimeV2;
    static double appendtime;
    static double SIMDparse1;
    static double SIMDparse2;
    static double SIMDparse3;
    static double SIMDcountbs;
    static double HYBDappendtime;
    bool googlevar = false;
};
#endif
#endif
