
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

#ifndef M3SYM_Database_h
#define M3SYM_Database_h

#include <unordered_set>

#include "common.h"

/// A collection class to hold instances of a given class

/*!
 *  The Database class is a template for holding a collection of objects
 *  It is used by all elements in the SubSystem to track instances of
 *  the given class. The databases are then used in various ways, i.e.
 *  mechanical minimization uses all beads for its Minimizer methods,
 *  ForceField uses the collections to calculate forces and energy, etc.
 !*/

template<class T>
class Database {
    
protected:
    unordered_set<T> _elems;  ///< The elements in the collection
    
    int _ID = 0; ///< Running unique index of each element
public:
    
    /// Add an element to the collection
    void addElement(T elem) {
        
        _elems.insert(elem);
    }
    
    /// Remove an element from the collection
    void removeElement(T elem) {
        
        //try to find
        auto it = _elems.find(elem);
        if(it != _elems.end()) _elems.erase(it);
    }
    
    /// Clear the contents of the database
    void clearElements() {
        
        _elems.clear();
    }
    
    /// Get all items in database
    const unordered_set<T>& getElements() {
        
        return _elems;
    }
    
    /// Count the number of objects in the collection
    int countElements() {
        
        return _elems.size();
    }
    
    ///Return a unique id
    int getID() { return _ID++;}
};


#endif
