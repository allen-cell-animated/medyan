//
//  ChemNRM.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemNRM_h
#define CytoSim_ChemNRM_h

#include "Reaction.h"

class ReactionNodeNRM {
private:
    Reaction *_react;
    float _tau;
    float _a;
public:
    ReactionNodeNRM(Reaction *RB) : _react (RB), _tau(std::numeric_limits<float>::quiet_NaN()) {
        std::cout << "ReactionNodeNRM ctor, ptr=" << this << std::endl;
        _a=0;
    }
    ReactionNodeNRM(const ReactionNodeNRM& other) = delete;
    ReactionNodeNRM& operator=(ReactionNodeNRM &rn) = delete;
    Reaction* getReaction() {return _react;};
    float getTau() const {return _tau;}
    void setTau(float tau){_tau=tau;}
    void makeStep(float t) {
        _react->makeStep();
        _a=_react->computePropensity();
        randDrawTau();
        updateNode();
        for(auto rit=_react->beginAffected();rit!=_react->endAffected();++rit){
            // get a_prev, tau_prev from the loop RNode
            // ...
            float a_prev=1.0; //fake
            float tau_prev=1.0;//fake
            float a_new = (*rit)->computePropensity();
            float tau_new = (a_prev/a_new)*(tau_prev-t)+t;
            // set tau_new ,a_new to the loop RNode
            // updateNode() the loop RNode
        }
    }
    void randDrawTau();
    void updateNode();
}; 

typedef ReactionNodeNRM RNode;

//
//template <typename PriorityQueue>
//class ChemNRM {
//    ChemNRM(){_init();}
//private:
//    void _init();
//    const PriorityQueue* _pq;
//};


#endif
