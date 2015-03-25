
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

#ifndef M3SYM_CMotorGhost_h
#define M3SYM_CMotorGhost_h

#include "common.h"

#include "CBound.h"

//FORWARD DECLARATIONS
class MotorGhost;
class SubSystem;
class CCylinder;
class Compartment;

/// A class to represent the chemical component of a MotorGhost.
/*!
 *  The CMotorGhost class contains chemical info of the parent MotorGhost.
 */

class CMotorGhost : public CBound {
    
private:
    MotorGhost* _pMotorGhost; ///< Pointer to parent
    
    CCylinder* _cc1; ///< Pointer to first CCylinder
    CCylinder* _cc2; ///< Pointer to second CCylinder
    
public:
    /// Default constructor and destructor
    CMotorGhost(Compartment* c, CCylinder* cc1, CCylinder* cc2)
        : CBound(c), _cc1(cc1), _cc2(cc2)  {}
    
    ~CMotorGhost() {}
    
    /// Copy constructor, standard
    CMotorGhost(const CMotorGhost& rhs, Compartment* c)
        : CBound(c), _pMotorGhost(rhs._pMotorGhost), _cc1(rhs._cc1), _cc2(rhs._cc2) {
        
        setFirstSpecies(rhs._firstSpecies);
        setSecondSpecies(rhs._secondSpecies);
        
        setOffReaction(rhs._offRxn);
    }
    
    /// Assignment is not allowed
    CMotorGhost& operator=(CMotorGhost &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CMotorGhost* clone(Compartment* c) {
        return new CMotorGhost(*this, c);
    }
    
    /// Set parent
    void setMotorGhost(MotorGhost* MotorGhost) {_pMotorGhost = MotorGhost;}
    /// Get parent
    MotorGhost* getMotorGhost() {return _pMotorGhost;}
    
    /// Create the off reaction for this MotorGhost
    void createOffReaction(ReactionBase* onRxn, float offRate, SubSystem* ps);
    
};

#endif
