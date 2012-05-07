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

class ReactionNodeNRM;
class PQNode;
typedef boost::heap::binomial_heap<PQNode> boost_heap;
typedef boost::heap::binomial_heap<PQNode>::handle_type handle_t;

class ReactionNodeNRM {
private:
    std::vector<ReactionNodeNRM*> _dependents;
    boost_heap& _heap;
    Reaction &_react;
    handle_t _handle;
public:
    ReactionNodeNRM(Reaction &r, boost_heap &bh);
    ~ReactionNodeNRM();
    ReactionNodeNRM(const ReactionNodeNRM& rhs) = delete;
    ReactionNodeNRM& operator=(ReactionNodeNRM &rhs) = delete;
    Reaction& getReaction() {return _react;};
    void updateHeap();
    float getTau() const;
    handle_t& getHandle();
    void setTau(float tau);
    void registerNewDependent(ReactionNodeNRM *rn); 
    void unregisterDependent(ReactionNodeNRM *rn);
    void makeStep(float t);
    void printSelf() const;
    void printDependents() const;
    //    void randDrawTau();
}; 


struct PQNode {
    ReactionNodeNRM *rnode;
    float tau;
    float a;
    bool operator<(PQNode const &rhs) const{
        return tau < rhs.tau;
    }
};

#endif
