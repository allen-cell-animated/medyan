
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

#ifndef M3SYM_GController_h
#define M3SYM_GController_h

#include <iostream>
#include <vector>

#include "common.h"

///OutOfBoundsException is a exception to be thrown when an index/coordinate is out of bounds of the grid
class OutOfBoundsException : public exception {
    
    virtual const char* what() const throw() {
        return "An element is out of the bounds of the grid.";
    }
};

///FORWARD DECLARATIONS
class Boundary;
class Compartment;

/// GController class is used to control the geometry of the grid, as well as the geometry of entire system
/*!
 *  The GeometryController class is used to control the geometry of the simulation, which includes the compartment 
 *  grid geometry as well as any boundary conditions that are in place. It has functionality to initialize a grid 
 *  based on the given dimensionality as well as find the correct compartment based on a set of coordinates.
 */
class GController {
    
private:
    ///local parameters stored for efficiency
    static short _nDim;
    static vector<int> _grid;
    static vector<double> _compartmentSize;
    
    ///Generate all neighbors lists for each compartment
    void generateConnections();
    
public:
    
    ///initialize the grid based on input parameters
    void initializeGrid();
    
    ///Activate compartments
    void activateCompartments(Boundary* boundary);
    
    /// Alternate getter from the grid
    static Compartment* getCompartment(const vector<size_t> &indices);
    /// Get the compartment given a set of coordinates
    static Compartment* getCompartment(const vector<double> &coords);
    
    /// Get all compartments within a given range from the specified coordinate
    /// @param ccheck - compartment to check when initially calling this function
    /// @param compartments - list of compartments that are within range. This will be populated by the function
    static void findCompartments(const vector<double>& coords, Compartment* ccheck, double dist, vector<Compartment*>& compartments);
    
};


#endif
