
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CCaMKIIingPoint_h
#define MEDYAN_CCaMKIIingPoint_h

#include "common.h"

#include "CBound.h"
#include "Compartment.h"

#define SPECIESB_BINDING_INDEX 0

//FORWARD DECLARATIONS
class CaMKIIingPoint;
class SubSystem;

/// A class to represent the chemical component of a CaMKIIingPoint.
/*!
 *  The CCaMKIIingPoint class contains chemical info of the parent CaMKIIingPoint.
 *
 *  Extending CBound, this class tracks its corresponding Species and unbinding Reaction.
 */
class CCaMKIIingPoint : public CBound {
    
private:
    CaMKIIingPoint* _pCaMKIIingPoint; ///< Pointer to parent
    
    short _camkiiType; ///< CaMKIIing point type
    
public:
    /// Default constructor and destructor
    /// @param pos - monomer index on first cylinder
    CCaMKIIingPoint(short camkiiType, Compartment* c,
                    CCylinder* cc1, CCylinder* cc2, int position);
    ///Destructor, removes off reaction from system
    ~CCaMKIIingPoint();
    
    /// Copy constructor, standard
    CCaMKIIingPoint(const CCaMKIIingPoint& rhs, Compartment* c)
        : CBound(rhs._filamentType, c, rhs._cc1, rhs._cc2, rhs._position1, rhs._position2),
          _pCaMKIIingPoint(rhs._pCaMKIIingPoint) {
        
          //set species
          setFirstSpecies(rhs._firstSpecies);
            
          //set reaction
          setOffReaction(rhs._offRxn);
            
          //set rates
          setOnRate(rhs._onRate);
          setOffRate(rhs._offRate);
    }
    
    /// Assignment is not allowed
    CCaMKIIingPoint& operator=(CCaMKIIingPoint &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CCaMKIIingPoint* clone(Compartment* c) {
        
        CCaMKIIingPoint* cb = new CCaMKIIingPoint(*this, c);
        _offRxn = nullptr; return cb;
    }
    
    /// Set parent
    void setCaMKIIingPoint(CaMKIIingPoint* CaMKIIingPoint) {
        _pCaMKIIingPoint = CaMKIIingPoint;
    }
    /// Get parent
    CaMKIIingPoint* getCaMKIIingPoint() {return _pCaMKIIingPoint;}
    
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps);
};

#endif
