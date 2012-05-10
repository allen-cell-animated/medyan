//
//  ChemNRMImpl.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/6/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemNRMImpl_h
#define CytoSim_ChemNRMImpl_h

#include <boost/heap/binomial_heap.hpp>

#include "Reaction.h"
#include "ChemNRM.h"

class PQNode;
class RNodeNRM;

typedef boost::heap::binomial_heap<PQNode> boost_heap;
typedef boost::heap::binomial_heap<PQNode>::handle_type handle_t;

class PQNode {
public: 
    RNodeNRM *_rn;
    double _tau;
    friend class RNodeNRM;
    PQNode(RNodeNRM *rnode) : _rn(rnode), _tau (std::numeric_limits<float>::quiet_NaN()) {}
    bool operator<(PQNode const &rhs) const{
        return _tau < rhs._tau;
    }
};

class RNode{
};

class RNodeNRM : public RNode {
private:
    handle_t _handle;
    Reaction *_react;
    float _a;
public:
    RNodeNRM(Reaction *r, boost_heap &heap);
    RNodeNRM(const RNodeNRM& rhs) = delete;
    RNodeNRM& operator=(RNodeNRM &rhs) = delete;
    Reaction* getReaction() const {return _react;};
    void updateHeap();
//    float getTau() const;
//    void setTau(float tau);
    float getTau() const {return (*_handle)._tau;}
    void setTau(float tau) {(*_handle)._tau=tau;}
    handle_t& getHandle() {return _handle;}
    void makeStep(float t);
    void printSelf() const;
    void printDependents() const;
    //    void randDrawTau();
};


#endif
