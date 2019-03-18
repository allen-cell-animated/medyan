
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

#include <cstddef> // size_t
#include <utility> // forward
#include <vector>

struct DatabaseDataDefault {
    void push_back() {}
    void pop_back() {}
    void copy_from_back(std::size_t dbIndex) {}
};

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
template< typename T, typename DatabaseData = DatabaseDataDefault >
class Database {
    
protected:
    static std::vector<T*> _elems;  ///< Pointer to the elements in the collection
    static DatabaseData _dbData;
    
    std::size_t _id;
    std::size_t _dbIndex;
    
    //DEPRECATED AS OF 9/22/16
//    int _transferID = -1; ///< index of a species ID to transfer
//                          ///< for now, this is used only in the case of motors.
//                          ///< If there is no transfer, the tag is marked as -1.

public:

    static const auto& getElements() { return _elems; }
    static auto& getDbData() { return _dbData; }

    // Get the next unique id
    static std::size_t nextId() { static std::size_t nextId = 0; return nextId++; }

    // Add element on construction
    template< typename... Args >
    Database(Args&&... args) : _id(nextId()), _dbIndex(_elems.size()) {
        _elems.push_back(static_cast<T*>(this));
        _dbData.push_back(std::forward<Args>(args)...);
    }
    // Remove element on destruction
    ~Database() {
        if(_dbIndex + 1 != _elems.size()) {
            // Move the data from the last element to the current position
            _elems[_dbIndex] = _elems.back();
            _dbData.copy_from_back(_dbIndex);
            // Updata _dbIndex of the original last element
            _elems[_dbIndex] -> _dbIndex = _dbIndex;
        }

        // Pop the last element
        _elems.pop_back();
        _dbData.pop_back();
    }

    std::size_t getId() const { return _id; }
    std::size_t getDbIndex() const { return _dbIndex; }
    
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

// Static variable definition (can be inlined starting C++17)
template< typename T, typename DatabaseData >
std::vector<T*> Database< T, DatabaseData >::_elems;
template< typename T, typename DatabaseData >
DatabaseData Database< T, DatabaseData>::_dbData;

template< typename T >
class OldDatabase {
protected:
    std::vector<T> _elems;  ///< Elements in the collection

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
    

};

#endif
