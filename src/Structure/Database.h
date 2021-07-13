
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
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

#include <algorithm> // max
#include <cstddef> // size_t
#include <utility> // forward
#include <vector>

#include <iostream>

//-----------------------------------------------------------------------------
// DatabaseData stores the data that has the same indexing with the internal
// Database structure.
// With stable indexing the following functions must be implemented:
//   void push_back(data...)
//   void set_content(size_t pos, data...)
//   void move_content(size_t from, size_t to)
//   void resize(size_t size)
// With dynamic indexing the following functions must be implemented:
//   void push_back(data...)
//   void pop_back()
//   void move_from_back(size_t pos)
//-----------------------------------------------------------------------------
struct DatabaseDataDefault {
    void push_back() {}
    void pop_back() {}
    void set_content(std::size_t) {}
    void move_from_back(std::size_t) {}
    void move_content(std::size_t from, std::size_t to) {}
    void resize(std::size_t) {}
    void settodummy(std::size_t) {}
};

template< typename T >
class DatabaseBase {
    
private:
    static std::vector<T*> _elems;  // Pointer to the elements in the collection
    static std::size_t _nextId;     // Next unique id

    std::size_t _id;
    std::size_t _index;

public:
    static const auto& getElements() { return _elems; }
    static auto numElements() { return _elems.size(); }

    // Add element on construction
    DatabaseBase() : _id(_nextId++), _index(_elems.size()) {
        _elems.push_back(static_cast<T*>(this));
    }
    // Remove element on destruction
    ~DatabaseBase() {
        if(_index + 1 != _elems.size()) {
            // Move the data from the last element to the current position
            _elems[_index] = _elems.back();
            // Updata _index of the original last element
            _elems[_index] -> DatabaseBase< T >::_index = _index;
        }

        // Pop the last element
        _elems.pop_back();
    }

    std::size_t getId() const { return _id; }
    std::size_t getIndex() const { return _index; }

    // This function overrides the current id of the element.
    // One should not use it unless in cases like re-initializing the system.
    void overrideId(std::size_t id) {
        _id = id;
        _nextId = std::max(_id + 1, _nextId);
    }

};
// Static variable definition (can be inlined starting C++17)
template< typename T >
std::vector<T*> DatabaseBase< T >::_elems;
template< typename T >
std::size_t DatabaseBase< T >::_nextId = 0;

template< typename DatabaseData > class DatabaseDataManager {

private:
    static DatabaseData _dbData;

public:
    static auto&       getDbData()      { return _dbData; }
    static const auto& getDbDataConst() { return _dbData; }
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

// Specialization for unstable indexing
template< typename T, typename DatabaseData >
class Database< T, false, DatabaseData > : public DatabaseBase<T>,
                                           public DatabaseDataManager<DatabaseData> {

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

// Specialization for stable indexing
template< typename T, typename DatabaseData >
class Database< T, true, DatabaseData > : public DatabaseBase<T>,
                                          public DatabaseDataManager<DatabaseData> {
    static std::vector<T*> _stableElems;
    static std::vector<std::size_t> _deletedIndices;

    std::size_t _stableIndex;

public:
    using db_base_type = DatabaseBase<T>;
    using db_data_type = DatabaseDataManager<DatabaseData>;

    static const auto& getStableElement(std::size_t pos) { return _stableElems[pos]; }
    // Get raw number of stable elements (including deleted)
    static auto rawNumStableElements() { return _stableElems.size(); }

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
        db_data_type::getDbData().settodummy(_stableIndex);
        _deletedIndices.push_back(_stableIndex);
    }

    std::size_t getStableIndex() const { return _stableIndex; }

    // Getting information for debug purposes
    static const std::vector<std::size_t>& getDeletedIndices() { return _deletedIndices; }

};
template< typename T, typename DatabaseData > std::vector<T*> Database< T, true, DatabaseData >::_stableElems;
template< typename T, typename DatabaseData > std::vector<std::size_t> Database< T, true, DatabaseData >::_deletedIndices;

#endif
