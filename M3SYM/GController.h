
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

#ifndef M3SYM_GController_h
#define M3SYM_GController_h

#include <vector>

#include "common.h"

/// An exception to be thrown when an index/coordinate is out of bounds of the grid
class OutOfBoundsException : public exception {
    
    virtual const char* what() const throw() {
        return "An element is out of the bounds of the grid. Try adjusting minimization parameters.";
    }
};

/// An exception to be thrown when an index/coordinate is NaN
class NaNCoordinateException : public exception {
    
    virtual const char* what() const throw() {
        return "A element coordinate is NaN. Try adjusting minimization parameters.";
    }
    
};

//FORWARD DECLARATIONS
class Boundary;
class Compartment;
class CompartmentGrid;

/// Used to control the geometry of the CompartmentGrid, as well as the geometry of
/// the entire system
/*!
 *  The GeometryController class is used by the SubSystem to control the geometry of 
 *  the simulation, which includes the CompartmentGrid geometry as well as any 
 *  [Boundaries] (@ref Boundary) that are in place. It has functionality to initialize
 *  a CompartmentGrid based on the given dimensionality as well as find the correct
 *  Compartment based on a set of coordinates.
 *  @note - The geometry controller currently only supports three-dimensional grids.
 */
class GController {
    
private:
    static short _nDim; ///< Number of dimensions in the system
    static vector<int> _grid; ///< Size of each dimension, in compartment lengths
    static vector<double> _compartmentSize; ///< Compartment size in each dimension
    static vector<double> _centerGrid; ///< The center of the grid
    
    static CompartmentGrid* _compartmentGrid; ///< The compartment grid
    
    ///Generate all neighbors lists for each compartment in the CompartmentGrid
    void generateConnections();
    
public:
    /// Initialize and return the grid based on input parameters
    CompartmentGrid* initializeGrid();
    
    /// Activate compartments in compartment grid based on a boundary
    void activateCompartments(Boundary* boundary);
    
    //@{
    /// Get a compartment based on coordinates or indices
    static Compartment* getCompartment(const vector<size_t> &indices);
    static Compartment* getCompartment(const vector<double> &coords);
    //@}
    
    /// Get the center of the grid space
    static const vector<double>& getCenter() {return _centerGrid;}

    /// Get all compartments within a given range from the specified coordinate
    /// @param ccheck - Compartment to check when initially calling this function
    /// @param compartments - List of compartments that are within range. This will be
    /// populated by the function
    static void findCompartments(const vector<double>& coords,
                                 Compartment* ccheck, double dist,
                                 vector<Compartment*>& compartments);
    
    /// Choose a random compartment from the grid (that is activated)
    static Compartment* getRandomCompartment();
    
    /// Get random coordinates in a given compartment
    static vector<double> getRandomCoordinates(Compartment* c);
};

#endif
