//
//  CLinker.h
//  Cyto
//
//  Created by James Komianos on 10/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CLinker__
#define __Cyto__CLinker__

#include <iostream>

#include "CBound.h"
#include "Compartment.h"

class Linker;

///CLinker is a class to represent the chemical component of a Linker
/*! 
 *  The CLinker class contains chemical info of the parent Linker.
 */

class CLinker : public CBound {
    
private:
    Linker* _pLinker; ///< ptr to parent linker

public:
    ///Default constructor and destructor
    CLinker(Compartment* c) :CBound(c, nullptr, nullptr) {}
    ~CLinker() {}
    
    ///Copy constructor, standard
    CLinker(const CLinker& rhs, Compartment* c) : _pLinker(rhs._pLinker), CBound(c, nullptr, nullptr) {}
    
    /// Assignment is not allowed
    CLinker& operator=(CLinker &rhs) = delete;
    
    ///Clone, calls copy constructor
    virtual CLinker* clone(Compartment* c) {
        return new CLinker(*this, c);
    }
    
    ///Setter and getter for parent linker
    void setLinker(Linker* linker) {_pLinker = linker;}
    Linker* getLinker() {return _pLinker;}
    
};


#endif /* defined(__Cyto__CLinker__) */
