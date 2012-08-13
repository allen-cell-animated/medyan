//
//  ChemRNode.h
//  CytoSim
//
//  Created by Garegin Papoian on 7/14/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemRNode_h
#define CytoSim_ChemRNode_h

namespace chem {

/// This is an abstract base class for classes that need to be associated with the given Reaction object.
class RNode{
public:
    /// Dtor
    virtual ~RNode() {}
    
    /// This method is called by Reaction::activateReaction(). Its effect depends on the underlying stochastic 
    /// simulatio algorithm. For example, in the NRM algorithm, a new tau and a are computed and the heap is updated. 
    virtual void activateReaction() = 0;

    /// This method is called by Reaction::passivateReaction(). Its effect depends on the underlying stochastic 
    /// simulatio algorithm. For example, in the NRM algorithm, a tau is set to infinity and the heap is updated. 
    virtual void passivateReaction() = 0;
    
    /// Return true if the Reaction is currently passivated
    virtual bool isPassivated() const = 0;
};

} // end of namespace


#endif
