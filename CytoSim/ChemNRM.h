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

class ReactionNode {
private:
    ReactionBase& _react;
    float _tau;
    const bool _forward;
public:
    ReactionNode(ReactionBase& RB, bool forward) : _react (RB), _tau(std::numeric_limits<float>::quiet_NaN()), _forward(forward) {};
    ReactionNode(const ReactionNode& other) : 
    _react (other._react), _tau(other._tau), _forward(other._forward) {};
    ReactionNode& operator=(ReactionNode &r) = delete;
    ReactionBase* getReactBase() {return &_react;};
    float getTau() const {return _tau;}
    void setTau(float tau){_tau=tau;}
    void makeStep() {
        if(_forward)
            _react.makeFStep();
        else
            _react.makeBStep();
        _react.updatePropensities();
    }
}; 



#endif
