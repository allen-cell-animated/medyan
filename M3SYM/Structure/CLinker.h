
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_CLinker_h
#define M3SYM_CLinker_h

#include "common.h"

#include "CBound.h"

//FORWARD DECLARATIONS
class Linker;
class SubSystem;
class CCylinder;
class Compartment;

/// To represent the chemical component of a Linker.
/*! 
 *  The CLinker class contains chemical info of the parent Linker.
 */

class CLinker : public CBound {
    
private:
    Linker* _pLinker; ///< Pointer to parent

public:
    ///Constructor
    ///@param pos1 - monomer index on first cylinder
    ///@param pos2 - monomer index on second cylinder
    CLinker(short linkerType, Compartment* c,
            CCylinder* cc1, CCylinder* cc2, int pos1, int pos2);
    
    ///Destructor, removes off reaction from system
    ~CLinker();
    
    /// Copy constructor, standard
    CLinker(const CLinker& rhs, Compartment* c)
        : CBound(c, rhs._cc1, rhs._cc2), _pLinker(rhs._pLinker){
        
        //set species
        setFirstSpecies(rhs._firstSpecies);
        setSecondSpecies(rhs._secondSpecies);
        
        //set reaction
        setOffReaction(rhs._offRxn);
            
        //set rates
        setOnRate(rhs._onRate);
        setOffRate(rhs._offRate);
    }
    
    /// Assignment is not allowed
    CLinker& operator=(CLinker &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CLinker* clone(Compartment* c) {
        CLinker* cl = new CLinker(*this, c);
        _offRxn = nullptr; return cl;
    }
    
    /// Set parent
    void setLinker(Linker* linker) {_pLinker = linker;}
    /// Get parent 
    Linker* getLinker() {return _pLinker;}
    
    /// Create the off reaction for this Linker
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps);
    
};


#endif