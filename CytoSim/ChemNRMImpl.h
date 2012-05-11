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
#include "ChemNRM.h"

class PQNode;
class RNodeNRM;
class ChemNRMImpl;

typedef boost::heap::pairing_heap<PQNode> boost_heap;
typedef boost::heap::pairing_heap<PQNode>::handle_type handle_t;

class PQNode {
public:
    PQNode(RNodeNRM *rnode) : _rn(rnode), _tau (std::numeric_limits<float>::quiet_NaN()) {}
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
};

class RNodeNRM : public RNode {
public:
    RNodeNRM(Reaction *r, boost_heap &heap);
    RNodeNRM(const RNodeNRM& rhs) = delete;
    RNodeNRM& operator=(RNodeNRM &rhs) = delete;
    Reaction* getReaction() const {return _react;};
    void updateHeap(boost_heap &heap);
    float getTau() const {return (*_handle)._tau;}
    void setTau(float tau) {(*_handle)._tau=tau;}
    handle_t& getHandle() {return _handle;}
    float getPropensity() const {return _a;}
    void reComputePropensity() {_a=_react->computePropensity ();}
    void makeStep() {_react->makeStep();}
    void printSelf() const;
    void printDependents() const;
private:
    handle_t _handle;
    Reaction *_react;
    float _a;
};

class ChemNRMImpl {
public:
    ChemNRMImpl() : 
    _eng(static_cast<unsigned long>(time(nullptr))), _exp_distr(0.0), _t(0.0), _n_reacts(0) {}
    ChemNRMImpl(const ChemNRMImpl &rhs) = delete;
    ChemNRMImpl& operator=(ChemNRMImpl &rhs) = delete;
    size_t getSize() const {return _n_reacts;}
    float getTime() const {return _t;}
    void addReaction(Reaction *r) {_map_rnodes.emplace(r,make_unique<RNodeNRM>(r,_heap)); ++_n_reacts;}
    void removeReaction(Reaction *r) {
        _map_rnodes.erase(r);
        --_n_reacts;
        //if we are here, then a request was for a reaction which was not found - a bug
        //throw std::runtime_error("Major bug: Trying to remove a reaction which is not registered.");
    }
    void initialize();
    void run(int steps) {
        for(int i=0; i<steps; ++i)
            _makeStep();
    }
    void printReactions() const;
private:
    void _makeStep();
    void _generateNewRandTau(RNodeNRM *rn);
private:
    std::unordered_map<Reaction*, std::unique_ptr<RNodeNRM>> _map_rnodes;
    boost_heap _heap;
    std::mt19937 _eng;
    std::exponential_distribution<float> _exp_distr;
    float _t;
    size_t _n_reacts;
};
#endif
