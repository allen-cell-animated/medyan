
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
 *  The GeometryController class is used to control the geometry of the simulation, which includes the [CompartmentGrid] 
 *  (@ref CompartmentGrid) geometry as well as any [Boundaries] (@ref Boundary) that are in place. It has functionality 
 *  to initialize a grid based on the given dimensionality as well as find the correct [Compartment] (@ref Compartment) 
 *  based on a set of coordinates.
 */
class GController {
    
private:
    static short _nDim; ///< Number of dimensions in the system
    static vector<int> _grid; ///< Size of each dimension, in compartment lengths
    static vector<double> _compartmentSize; ///< Compartment size in each dimension
    
    ///Generate all neighbors lists for each compartment in the [CompartmentGrid] (@ref CompartmentGrid)
    void generateConnections();
    
public:
    
    /// Initialize the grid based on input parameters
    void initializeGrid();
    
    /// Activate compartments based on a [Boundary] (@ref Boundary)
    void activateCompartments(Boundary* boundary);
    
    //@{
    /// Get a [Compartment] (@ref Compartment) based on coordinates or indices
    static Compartment* getCompartment(const vector<size_t> &indices);
    static Compartment* getCompartment(const vector<double> &coords);
    //@}

    /// Get all [Compartments] (@ref Compartment) within a given range from the specified coordinate
    /// @param ccheck - Compartment to check when initially calling this function
    /// @param compartments - List of compartments that are within range. This will be populated by the function
    static void findCompartments(const vector<double>& coords, Compartment* ccheck,
                                 double dist, vector<Compartment*>& compartments);
    
};

#endif
