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

//class ReactionNode {
//private:
//    ReactionBase& _react;
//    float _tau;
//    const bool _forward;
//public:
//    ReactionNode(ReactionBase& RB, bool forward) : _react (RB), _tau(std::numeric_limits<float>::quiet_NaN()), _forward(forward) {
//        std::cout << "ReactionNode ctor, ptr=" << this << std::endl;
//    }
//    ReactionNode(const ReactionNode& other) : 
//    _react (other._react), _tau(other._tau), _forward(other._forward) {
//        std::cout << "ReactionNode copy ctor, ptr=" << this << std::endl;
//    };
//    ReactionNode& operator=(ReactionNode &rn) = delete;
//    ReactionBase* getReactBase() {return &_react;};
//    float getTau() const {return _tau;}
//    void setTau(float tau){_tau=tau;}
//    void makeStep() {
//        if(_forward)
//            _react.makeFStep();
//        else
//            _react.makeBStep();
//        _react.updatePropensities();
//    }
//}; 
//
//template <typename PriorityQueue>
//class ChemNRM {
//    ChemNRM(){_init();}
//private:
//    void _init();
//    const PriorityQueue* _pq;
//};


#endif
