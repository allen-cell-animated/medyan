
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_NeighborListContainer_h
#define M3SYM_NeighborListContainer_h

#include "common.h"

#include "NeighborListDB.h"

/// CylinderNLContainer class holds a [CylinderNeighborList](@ref CylinderNeighborList).
class CylinderNLContainer {

private:
    CylinderNeighborList* _neighborList;
    
public:
    /// Constructor, adds a cylinder neighbor list to the database
    CylinderNLContainer(float rMax = 0.0, float rMin = 0.0, bool crossFilamentOnly = false) {
        
        _neighborList = NeighborListDB::instance()->
                        createCylinderNeighborList(rMax, rMin, crossFilamentOnly);
    }
    /// Destructor, removes cylinder neighbor list from the database
    ~CylinderNLContainer() { NeighborListDB::instance()->removeNeighborList(_neighborList); }
    
    /// Get neighbor list
    CylinderNeighborList* getNeighborList() {return _neighborList;}
    
};

/// BoundaryElementNLContainer class holds a [BoundaryElementNeighborList](@ref BoundaryElementNeighborList).
class BoundaryElementNLContainer {
    
private:
    BoundaryElementNeighborList* _neighborList;

public:
    /// Constructor, adds a cylinder neighbor list to the database
    BoundaryElementNLContainer(float rMax = 0.0, float rMin = 0.0) {
        
        _neighborList = NeighborListDB::instance()->
                        createBoundaryElementNeighborList(rMax, rMin);
    }
    /// Destructor, removes boundary element neighbor list from the database
    ~BoundaryElementNLContainer() { NeighborListDB::instance()->removeNeighborList(_neighborList); }
    
    /// Get neighbor list
    BoundaryElementNeighborList* getNeighborList() {return _neighborList;}
};


#endif
