
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
 *  The Database also contains a holder for a transfer ID of any species.
 *  This is used when diffusing species ID's must be tracked, and these
 *  ID's are accessed by the respective binding managers.
 *
 *  @param T - class to hold
 *
 */
template<class T>
class Database {
    
protected:
    vector<T> _elems;  ///< The elements in the collection
    
    int _ID = 0; ///< Running unique index of each element
    
    //DEPRECATED AS OF 9/22/16
//    int _transferID = -1; ///< index of a species ID to transfer
//                          ///< for now, this is used only in the case of motors.
//                          ///< If there is no transfer, the tag is marked as -1.
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
    
    ///Used for a deletion of ID
    int deleteID() {return --_ID;}
    
    //DEPRECATED AS OF 9/22/16
//
//    //@{
//    ///Setters and getters for transfer ID
//    void setTransferID(int ID) {_transferID = ID;}
//    
//    int getTransferID() {
//        
//        int retID = _transferID;
//        _transferID = -1;
//        return retID;
//    }
    
    //@}
};


#endif
