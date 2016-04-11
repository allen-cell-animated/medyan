
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Database_h
#define MEDYAN_Database_h

#include <set>

#include "common.h"

/// A collection class to hold instances of a given class

/*!
 *  The Database class is a template for holding a collection of objects
 *  It is used by all elements in the SubSystem to track instances of
 *  the given class. The databases are then used in various ways, i.e.
 *  mechanical minimization uses all beads for its Minimizer methods,
 *  ForceField uses the collections to calculate forces and energy, etc.
 *  
 *  @param T - class to hold
 *
 */
template<class T>
class Database {
    
protected:
    vector<T> _elems;  ///< The elements in the collection
    
    int _ID = 0; ///< Running unique index of each element
public:
    
    /// Add an element to the collection
    void addElement(T elem) {
        
        _elems.push_back(elem);
    }
    
    /// Remove an element from the collection
    void removeElement(T elem) {
        
        //try to find
        auto it = find(_elems.begin(), _elems.end(), elem);
        if(it != _elems.end()) _elems.erase(it);
    }
    
    /// Clear the contents of the database
    void clearElements() { _elems.clear(); }
    
    /// Get all items in database
    vector<T>& getElements() { return _elems; }
    
    /// Count the number of objects in the collection
    int countElements() { return _elems.size(); }
    
    ///Return a unique id
    int getID() { return _ID++;}
};


#endif
