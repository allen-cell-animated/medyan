
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

#ifndef M3SYM_CLinker_h
#define M3SYM_CLinker_h

#include "common.h"

#include "CBound.h"

//FORWARD DECLARATIONS
class Linker;
class Compartment;

/// To represent the chemical component of a [Linker](@ref Linker).
/*! 
 *  The CLinker class contains chemical info of the parent [Linker](@ref Linker).
 */

class CLinker : public CBound {
    
private:
    Linker* _pLinker; ///< Pointer to parent linker

public:
    CLinker(Compartment* c) :CBound(c) {}
    ~CLinker() {}
    
    /// Copy constructor, standard
    CLinker(const CLinker& rhs, Compartment* c) : _pLinker(rhs._pLinker), CBound(c) {
        
        setFirstSpecies(rhs._firstSpecies);
        setSecondSpecies(rhs._secondSpecies);
    }
    
    /// Assignment is not allowed
    CLinker& operator=(CLinker &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CLinker* clone(Compartment* c) {
        return new CLinker(*this, c);
    }
    
    /// Set parent linker
    void setLinker(Linker* linker) {_pLinker = linker;}
    /// Get parent linker
    Linker* getLinker() {return _pLinker;}
    
};


#endif