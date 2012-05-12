//
//  ChemNRMImpl.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/6/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemNRMImpl_h
#define CytoSim_ChemNRMImpl_h

#include <vector>
#include <random>

#include <boost/heap/pairing_heap.hpp>

#include "utility.h"
#include "Reaction.h"
#include "ChemSim.h"

class PQNode;
class RNodeNRM;
class ChemNRMImpl;

typedef boost::heap::pairing_heap<PQNode> boost_heap;
typedef boost::heap::pairing_heap<PQNode>::handle_type handle_t;

class PQNode {
public:
    PQNode(RNodeNRM *rnode) : _rn(rnode), _tau (std::numeric_limits<double>::infinity()) {}
    ~PQNode(){_rn=nullptr;}
    bool operator<(PQNode const &rhs) const{
        return _tau > rhs._tau;
    }
private: 
    RNodeNRM *_rn;
    double _tau;
    friend class ChemNRMImpl;
    friend class RNodeNRM;
};

class RNode{
public:
    virtual ~RNode() {}
    virtual void activateReaction() = 0;
    virtual void passivateReaction() = 0;    
};

class RNodeNRM : public RNode {
public:
    RNodeNRM(Reaction *r, ChemNRMImpl &chem_nrm);
    RNodeNRM(const RNodeNRM& rhs) = delete;
    RNodeNRM& operator=(RNodeNRM &rhs) = delete;
    ~RNodeNRM(); 
    void generateNewRandTau();
    Reaction* getReaction() const {return _react;};
    void updateHeap();
    double getTau() const {return (*_handle)._tau;}
    void setTau(double tau) {(*_handle)._tau=tau;}
    handle_t& getHandle() {return _handle;}
    double getPropensity() const {return _a;}
    int getReactantsProduct() {return _react->getReactantsProduct();};
    void reComputePropensity() {_a=_react->computePropensity ();}
    void makeStep() {_react->makeStep();}
    void activateReaction();
    void passivateReaction();   
    void printSelf() const;
    void printDependents() const;
private:
    ChemNRMImpl &_chem_nrm;
    handle_t _handle;
    Reaction *_react;
    double _a;
};

class ChemNRMImpl : public ChemSimImpl {
public:
    ChemNRMImpl() : 
    ChemSimImpl(), _eng(static_cast<unsigned long>(time(nullptr))), _exp_distr(0.0), _t(0.0), _n_reacts(0) {}
    ChemNRMImpl(const ChemNRMImpl &rhs) = delete;
    ChemNRMImpl& operator=(ChemNRMImpl &rhs) = delete;
    ~ChemNRMImpl();
    size_t getSize() const {return _n_reacts;}
    double getTime() const {return _t;}
    boost_heap* getHeap() {return &_heap;} 
    void addReaction(Reaction *r);
    void removeReaction(Reaction *r);
    double generateTau(double a);
    void initialize();
    void run(int steps) {
        for(int i=0; i<steps; ++i){
            makeStep();
            if(i%1000000==0)
                std::cout << "ChemNRMImpl::run(): i=" << i << std::endl;
        }
    }
    void printReactions() const;
private:
    void makeStep();
private:
    std::unordered_map<Reaction*, std::unique_ptr<RNodeNRM>> _map_rnodes;
    boost_heap _heap;
    std::mt19937 _eng;
    std::exponential_distribution<double> _exp_distr;
    double _t;
    size_t _n_reacts;
};
#endif
