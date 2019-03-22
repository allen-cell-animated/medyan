
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
    void copy_from_back(std::size_t index) {}
    void move_content(std::size_t from, std::size_t to) {}
    void resize(std::size_t) {}
};

template< typename T >
class DatabaseBase {
    
private:
    static std::vector<T*> _elems;  ///< Pointer to the elements in the collection

    std::size_t _id;
    std::size_t _index;

public:
    static const auto& getElements() { return _elems; }
    static auto numElements() { return _elems.size(); }

    // Get the next unique id
    static std::size_t nextId() { static std::size_t nextId = 0; return nextId++; }

    // Add element on construction
    DatabaseBase() : _id(nextId()), _index(_elems.size()) {
        _elems.push_back(static_cast<T*>(this));
    }
    // Remove element on destruction
    ~DatabaseBase() {
        if(_index + 1 != _elems.size()) {
            // Move the data from the last element to the current position
            _elems[_index] = _elems.back();
            // Updata _index of the original last element
            _elems[_index] -> _index = _index;
        }

        // Pop the last element
        _elems.pop_back();
    }

    std::size_t getId() const { return _id; }
    std::size_t getIndex() const { return _index; }

};
// Static variable definition (can be inlined starting C++17)
template< typename T >
std::vector<T*> DatabaseBase< T >::_elems;

template< typename DatabaseData > class DatabaseDataManager {

private:
    static DatabaseData _dbData;

public:
    static auto& getDbData() { return _dbData; }
};
// Static variable definition
template< typename DatabaseData > DatabaseData DatabaseDataManager< DatabaseData >::_dbData;


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
template< typename T, bool stableIndexing, typename DatabaseData = DatabaseDataDefault > class Database;
template< typename T, typename DatabaseData >
class Database< T, false, DatabaseData > : public DatabaseBase<T>, public DatabaseDataManager<DatabaseData> {

public:

    // Add element on construction
    template< typename... Args >
    Database(Args&&... args) {
        getDbData().push_back(std::forward<Args>(args)...);
    }
    // Remove element on destruction (by swapping)
    ~Database() {
        if(getIndex() + 1 != getElements().size()) {
            // Move the data from the last element to the current position
            getDbData().copy_from_back(getIndex());
        }

        // Pop the last element
        getDbData().pop_back();
    }
};
template< typename T, typename DatabaseData >
class Database< T, true, DatabaseData > : public DatabaseBase<T>, public DatabaseDataManager<DatabaseData> {
    static std::vector<T*> _stableElems;
    static std::vector<std::size_t> _deletedIndices;

    std::size_t _stableIndex;

public:
    static const auto& getStableElement(std::size_t pos) { return _stableElems[pos]; }

    static void rearrange() {
        using std::size_t;

        const size_t numDeleted = _deletedIndices.size();
        const size_t currentSize = _stableElems.size();
        const size_t finalSize = currentSize - numDeleted;

        std::vector<char> isDeleted(numDeleted);
        // Mark to-be-deleted items with indices bigger than finalSize as deleted
        for(size_t i = 0; i < numDeleted; ++i)
            if(_deletedIndices[i] >= finalSize)
                isDeleted[_deletedIndices[i] - finalSize] = true;

        // Move the not-to-be-deleted items with bigger indices to the to-be-deleted items with small indices
        for(size_t indAfterFinal = 0, i = 0; indAfterFinal < numDeleted; ++indAfterFinal) {
            if(!isDeleted[indAfterFinal]) {
                while(i < numDeleted && _deletedIndices[i] >= finalSize) ++i; // Find (including current i) the next i with small index
                if(i < numDeleted) {
                    // Found. This should always be satisfied.
                    _stableElems[_deletedIndices[i]] = _stableElems[finalSize + indAfterFinal];
                    getDbData().move_content(finalSize + indAfterFinal, _deletedIndices[i]);
                    _stableElems[_deletedIndices[i]]->_stableIndex = _deletedIndices[i];
                }
            }
        }

        // Remove garbage
        _stableElems.resize(finalSize);
        getDbData().resize(finalSize);
        _deletedIndices.clear();
    }

    template< typename... Args >
    Database(Args&&... args) : _stableIndex(_stableElems.size()) {
        getDbData().push_back(std::forward<Args>(args)...);
    }
    ~Database() {
        // Only mark as deleted
        _deletedIndices.push_back(_stableIndex);
    }

    std::size_t getStableIndex() const { return _stableIndex; }
};
template< typename T, typename DatabaseData > std::vector<T*> Database< T, true, DatabaseData >::_stableElems;
template< typename T, typename DatabaseData > std::vector<std::size_t> Database< T, true, DatabaseData >::_deletedIndices;


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
    std::vector<T>& getElements() { return _elems; }
    
    /// Count the number of objects in the collection
    int countElements() { return _elems.size(); }
    

};

#endif
