//
//  CMotorGhost.h
//  Cyto
//
//  Created by James Komianos on 10/20/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CMotorGhost__
#define __Cyto__CMotorGhost__

#include <iostream>
#include "CBound.h"
#include "Compartment.h"

class MotorGhost;

///CMotorGhost is a class to represent the chemical component of a motor
/*!
 *  The CMotorGhost class contains chemical info of the parent MotorGhost.
 */

class CMotorGhost : public CBound {
    
private:
    MotorGhost* _pMotorGhost; ///< ptr to parent motorghost
    
public:
    ///Default constructor and destructor
    CMotorGhost(Compartment* c) :CBound(c) {}
    ~CMotorGhost() {}
    
    ///Copy constructor, standard
    CMotorGhost(const CMotorGhost& rhs, Compartment* c) : _pMotorGhost(rhs._pMotorGhost), CBound(c) {
        
        setFirstSpecies(rhs._firstSpecies);
        setSecondSpecies(rhs._secondSpecies);
    }
    
    /// Assignment is not allowed
    CMotorGhost& operator=(CMotorGhost &rhs) = delete;
    
    ///Clone, calls copy constructor
    virtual CMotorGhost* clone(Compartment* c) {
        return new CMotorGhost(*this, c);
    }
    
    ///Setter and getter for parent linker
    void setMotorGhost(MotorGhost* MotorGhost) {_pMotorGhost = MotorGhost;}
    MotorGhost* getMotorGhost() {return _pMotorGhost;}
    
};



#endif /* defined(__Cyto__CMotorGhost__) */
