
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

#ifndef M3SYM_NeighborListContainer_h
#define M3SYM_NeighborListContainer_h

#include "common.h"

#include "NeighborListImpl.h"
#include "NeighborListDB.h"

/// Holds a CCNeighborList.
class CCNLContainer {

protected:
    CCNeighborList* _neighborList;
    
public:
    /// Constructor, adds a cylinder-cylinder neighbor list to the database
    CCNLContainer(float rMax = 0.0, float rMin = 0.0) {
        
        _neighborList = new CCNeighborList(rMax, rMin);
        NeighborListDB::instance()->addNeighborList(_neighborList);
    }
    /// Destructor, removes cylinder-cylinder neighbor list from the database
    ~CCNLContainer() {
        NeighborListDB::instance()->removeNeighborList(_neighborList);
        delete _neighborList;
    }
    
    /// Get neighbor list
    CCNeighborList* getNeighborList() {return _neighborList;}
    
};

/// Holds a BBENeighborList.
class BBENLContainer {
    
protected:
    BBENeighborList* _neighborList;

public:
    /// Constructor, adds a bead-boundary element neighbor list to the database
    BBENLContainer(float rMax = 0.0) {
        
        _neighborList = new BBENeighborList(rMax);
        NeighborListDB::instance()->addNeighborList(_neighborList);
    }
    /// Destructor, removes bead - boundary element neighbor list from the database
    ~BBENLContainer() {
        NeighborListDB::instance()->removeNeighborList(_neighborList);
        delete _neighborList;
    }
    
    /// Get neighbor list
    BBENeighborList* getNeighborList() {return _neighborList;}
};

/// Holds a CBENeighborList.
class CBENLContainer {
    
protected:
    CBENeighborList* _neighborList;
    
public:
    /// Constructor, adds a cylinder - boundary element neighbor list to the database
    CBENLContainer(float rMax = 0.0) {
        
        _neighborList = new CBENeighborList(rMax);
        NeighborListDB::instance()->addNeighborList(_neighborList);
    }
    /// Destructor, removes cylinder - boundary element neighbor list from the database
    ~CBENLContainer() {
        NeighborListDB::instance()->removeNeighborList(_neighborList);
        delete _neighborList;
    }
    
    /// Get neighbor list
    CBENeighborList* getNeighborList() {return _neighborList;}
};


#endif
