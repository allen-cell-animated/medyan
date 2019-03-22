
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

//-----------------------------------------------------------------------------
// DatabaseData stores the data that has the same indexing with the internal
// Database structure.
// With stable indexing the following functions must be implemented:
//   push_back(data...)
//   set_content(size_t pos, data...)
//   move_content(size_t from, size_t to)
//   resize(size_t size)
// With dynamic indexing the following functions must be implemented:
//   push_back(data...)
//   pop_back()
//   move_from_back(size_t pos)
//-----------------------------------------------------------------------------
struct DatabaseDataDefault {
    void push_back() {}
    void pop_back() {}
    void set_content(std::size_t) {}
    void move_from_back(std::size_t) {}
    void move_content(std::size_t from, std::size_t to) {}
    void resize(std::size_t) {}
};

template< typename T >
class DatabaseBase {
    
private:
    static std::vector<T*> _elems;  ///< Pointer to the elements in the collection

    std::size_t _id;
    std::size_t _index;

    // Get the next unique id
    static std::size_t _nextId() { static std::size_t nextId = 0; return nextId++; }

public:
    static const auto& getElements() { return _elems; }
    static auto numElements() { return _elems.size(); }

    // Add element on construction
    DatabaseBase() : _id(_nextId()), _index(_elems.size()) {
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
 *  @param T - class to point to
 *  @param stableIndexing - when this is set to true, the database will be able
 *    to track an additional indexing which will never be invalidated until
 *    rearrange() is called.
 *  @param DatabaseData - the managed vectorized data. When stableIndexing is
 *    true, one uses getStableIndex() to access the data, but the size of the
 *    data might be bigger than the number of elements. When stableIndexing is
 *    false, one uses getIndex() to access the data, where the index might be
 *    invalidated after any other element is created/destroyed.
 */
template< typename T, bool stableIndexing, typename DatabaseData = DatabaseDataDefault > class Database;
template< typename T, typename DatabaseData >
class Database< T, false, DatabaseData > : public DatabaseBase<T>, public DatabaseDataManager<DatabaseData> {

public:

    using db_base_type = DatabaseBase<T>;
    using db_data_type = DatabaseDataManager<DatabaseData>;

    // Add element on construction
    template< typename... Args >
    Database(Args&&... args) {
        db_data_type::getDbData().push_back(std::forward<Args>(args)...);
    }
    // Remove element on destruction (by swapping)
    ~Database() {
        if(this->getIndex() + 1 != db_base_type::getElements().size()) {
            // Move the data from the last element to the current position
            db_data_type::getDbData().move_from_back(this->getIndex());
        }

        // Pop the last element
        db_data_type::getDbData().pop_back();
    }
};
template< typename T, typename DatabaseData >
class Database< T, true, DatabaseData > : public DatabaseBase<T>, public DatabaseDataManager<DatabaseData> {
    static std::vector<T*> _stableElems;
    static std::vector<std::size_t> _deletedIndices;

    std::size_t _stableIndex;

public:
    using db_base_type = DatabaseBase<T>;
    using db_data_type = DatabaseDataManager<DatabaseData>;

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
                    db_data_type::getDbData().move_content(finalSize + indAfterFinal, _deletedIndices[i]);
                    _stableElems[_deletedIndices[i]]->_stableIndex = _deletedIndices[i];
                }
                ++i;
            }
        }

        // Remove garbage
        _stableElems.resize(finalSize);
        db_data_type::getDbData().resize(finalSize);
        _deletedIndices.clear();
    }

    template< typename... Args >
    Database(Args&&... args) {
        if(_deletedIndices.empty()) {
            // Append at the back
            _stableIndex = _stableElems.size();
            _stableElems.push_back(static_cast<T*>(this));
            db_data_type::getDbData().push_back(std::forward<Args>(args)...);
        } else {
            // Fill in the last hole
            _stableIndex = _deletedIndices.back();
            _stableElems[_stableIndex] = static_cast<T*>(this);
            db_data_type::getDbData().set_content(_stableIndex, std::forward<Args>(args)...);
            _deletedIndices.pop_back();
        }
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
