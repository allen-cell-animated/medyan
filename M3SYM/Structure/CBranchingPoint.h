
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

#ifndef M3SYM_CBranchingPoint_h
#define M3SYM_CBranchingPoint_h

#include "common.h"

#include "CBound.h"
#include "Compartment.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// A class to represent the chemical component of a BranchingPoint.
/*!
 *  The CBranchingPoint class contains chemical info of the parent BranchingPoint.
 */
class CBranchingPoint : public CBound {
    
private:
    BranchingPoint* _pBranchingPoint; ///< Pointer to parent
    
public:
    /// Default constructor and destructor
    CBranchingPoint(Compartment* c) :CBound(c) {}
    ~CBranchingPoint() {}
    
    /// Copy constructor, standard
    CBranchingPoint(const CBranchingPoint& rhs, Compartment* c)
        : _pBranchingPoint(rhs._pBranchingPoint), CBound(c) {
        
          setFirstSpecies(rhs._firstSpecies);
          //setOffReaction(rhs._offRxn);
    }
    
    /// Assignment is not allowed
    CBranchingPoint& operator=(CBranchingPoint &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CBranchingPoint* clone(Compartment* c) {
        return new CBranchingPoint(*this, c);
    }
    
    /// Set parent
    void setBranchingPoint(BranchingPoint* BranchingPoint)
        {_pBranchingPoint = BranchingPoint;}
    /// Get parent
    BranchingPoint* getBranchingPoint() {return _pBranchingPoint;}
};

#endif
