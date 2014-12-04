
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

#ifndef M3SYM_CBranchPoint_h
#define M3SYM_CBranchPoint_h

#include "common.h"

#include "CBound.h"
#include "Compartment.h"

//FORWARD DECLARATIONS
class BranchPoint;

/// A class to represent the chemical component of a BranchPoint.
/*!
 *  The CBranchPoint class contains chemical info of the parent BranchPoint.
 */

class CBranchPoint : public CBound {
    
private:
    BranchPoint* _pBranchPoint; ///< Pointer to parent
    
public:
    /// Default constructor and destructor
    CBranchPoint(Compartment* c) :CBound(c) {}
    ~CBranchPoint() {}
    
    /// Copy constructor, standard
    CBranchPoint(const CBranchPoint& rhs, Compartment* c) : _pBranchPoint(rhs._pBranchPoint), CBound(c) {
        
//        setFirstSpecies(rhs._firstSpecies);
//        setSecondSpecies(rhs._secondSpecies);
    }
    
    /// Assignment is not allowed
    CBranchPoint& operator=(CBranchPoint &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CBranchPoint* clone(Compartment* c) {
        return new CBranchPoint(*this, c);
    }
    
    /// Set parent
    void setBranchPoint(BranchPoint* branchPoint) {_pBranchPoint = branchPoint;}
    /// Get parent
    BranchPoint* getBranchPoint() {return _pBranchPoint;}
};

#endif
