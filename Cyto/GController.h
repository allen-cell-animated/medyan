//
//  GController.h
//  Cyto
//
//  Created by James Komianos on 8/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__GController__
#define __Cyto__GController__

#include <iostream>
#include <vector>
#include "CompartmentContainer.h"
#include "SystemParameters.h"

///Exception to be thrown when an index/coordinate is out of bounds of the grid
class OutOfBoundsException : public std::exception {
    
    virtual const char* what() const throw() {
        return "An element is out of the bounds of the grid.";
    }
};

class Boundary;

/// GController class is used to control the geometry of the grid, as well as the geometry of entire system
/*!
 *
 *  The GeometryController class is used by the SubSystem class to control the geometry of the simulation.
 *  This includes the compartment grid geometry as well as any boundary conditions that are in place.
 *  It has functionality to initialize a grid based on the given dimensionality as well as find the correct
 *  compartment based on a set of coordinates.
 */

class GController {
    
private:
    ///local parameters stored for efficiency
    static short _nDim;
    static std::vector<int> _grid;
    static std::vector<double> _compartmentSize;
    
    ///Generate all neighbors lists for each compartment
    void generateConnections();
    
public:
    
    ///initialize the grid based on input parameters
    void initializeGrid();
    
    ///Activate compartments
    void activateCompartments(Boundary* boundary);
    
    /// Alternate getter from the grid
    static Compartment* getCompartment(const std::vector<size_t> &indices)
    {
        size_t index = 0;
        size_t i = 0;
        for(auto x: indices)
        {
            
            ///Flatten the indices to 1D
            if(i == 0) {
                if(x >= _grid[0]) { throw OutOfBoundsException();}
                index += x;
            }
            else if(i == 1) {
                if(x >= _grid[1]) { throw OutOfBoundsException();}
                index += x * _grid[1];
            }
            else {
                if(x >= _grid[2]) { throw OutOfBoundsException();}
                index += x * _grid[0] * _grid[1];
            }
            
            i++;
        }
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
    /// Get the compartment given a set of coordinates
    static Compartment* getCompartment(const std::vector<double> &coords)
    {
        ///Check if out of bounds
        size_t index = 0;
        size_t i = 0;
        for(auto x: coords)
        {
            ///Flatten the coordinates to 1D, get integer index
            if(i == 0) {
                if(x < 0 || x >= (_compartmentSize[0] * _grid[0])) {
                    throw OutOfBoundsException();
                }
                index += int(x / _compartmentSize[0]);
            }
            else if(i == 1) {
                if(x < 0 || x >= (_compartmentSize[1] * _grid[1])) {
                    throw OutOfBoundsException();
                }
                index += int(x / _compartmentSize[1]) * _grid[0];
            }
            else {
                if(x < 0 || x >= (_compartmentSize[2] * _grid[2])) {
                    throw OutOfBoundsException();
                }
                index += int(x / _compartmentSize[2]) * _grid[0] * _grid[1];
            }
            i++;
        }
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
    /// Get all compartments within a given range from the specified compartment
    /// @param ccheck - compartment to check. when initially calling this function, ccheck should be the same as c
    /// @param compartments - list of compartments that are within range. This will be populated by the function
    static void findCompartments(const std::vector<double>& coords, Compartment* ccheck, double dist, std::vector<Compartment*>& compartments);
    
};


#endif /* defined(__Cyto__GController__) */
