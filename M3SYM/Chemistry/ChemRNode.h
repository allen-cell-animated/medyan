
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_ChemRNode_h
#define M3SYM_ChemRNode_h

#include "common.h"

/// This is an abstract base class for classes that need to be associated with the given [Reaction] (@ref Reaction) object.
class RNode{
public:
    /// Dtor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~RNode() noexcept {}
    
    /// This method is called by Reaction::activateReaction(). Its effect depends on the underlying stochastic 
    /// simulatio algorithm. For example, in the NRM algorithm, a new tau and a are computed and the heap is updated. 
    virtual void activateReaction() = 0;

    /// This method is called by Reaction::passivateReaction(). Its effect depends on the underlying stochastic 
    /// simulatio algorithm. For example, in the NRM algorithm, a tau is set to infinity and the heap is updated. 
    virtual void passivateReaction() = 0;
    
    /// Return true if the Reaction is currently passivated
    virtual bool isPassivated() const = 0;
};

#endif