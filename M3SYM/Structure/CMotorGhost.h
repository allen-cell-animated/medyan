
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

#ifndef M3SYM_CMotorGhost_h
#define M3SYM_CMotorGhost_h

#include "common.h"

#include "CBound.h"
#include "Compartment.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// A class to represent the chemical component of a MotorGhost.
/*!
 *  The CMotorGhost class contains chemical info of the parent MotorGhost.
 */

class CMotorGhost : public CBound {
    
private:
    MotorGhost* _pMotorGhost; ///< Pointer to parent
    
public:
    /// Default constructor and destructor
    CMotorGhost(Compartment* c) :CBound(c) {}
    ~CMotorGhost() {}
    
    /// Copy constructor, standard
    CMotorGhost(const CMotorGhost& rhs, Compartment* c) : _pMotorGhost(rhs._pMotorGhost), CBound(c) {
        
        setFirstSpecies(rhs._firstSpecies);
        setSecondSpecies(rhs._secondSpecies);
        
        _offRxn = rhs._offRxn;
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
    
};

#endif
