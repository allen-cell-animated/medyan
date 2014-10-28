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
    static short _nDim; ///<grid dimensionality
    static std::vector<int> _grid; ///< grid dimensions (in units of compartments)
    static std::vector<double> _compartmentSize; ///< compartment size in nm
    
    ///Generate all neighbors lists for each compartment
    static void generateConnections();
    
public:
    
    ///initialize the grid based on input parameters
    static void initializeGrid();
    
    ///Activate compartments
    static void activateCompartments(Boundary* boundary);
    
    /// Get compartment from the grid
    /// @param - args, the indices in n-dimensions of the compartment
//    template<typename ...Args>
//    static Compartment* getCompartment(Args&& ...args)
//    {
//        size_t index = 0;
//        size_t i = _nDim-1;
//        for(auto x: {args...})
//        {
//            index+= x * std::pow(_grid[i],i);
//            --i;
//        }
//        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
//        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
//    }
    
    /// Alternate getter from the grid
    static Compartment* getCompartment(const std::vector<size_t> &indices)
    {
        size_t index = 0;
        size_t i = 0;
        for(auto x: indices)
        {
            if(x >= _grid[i]) { throw OutOfBoundsException();}
            
            ///Flatten the indices to 1D
            if(i == 0)
                index += x;
            else if(i == 1)
                index += x * _grid[0];
            else
                index += x * _grid[0] * _grid[1];
            
            i++;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
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
            if(x < 0 || x >= (_compartmentSize[i] * _grid[i])) {
                throw OutOfBoundsException();
            }
            
            ///Flatten the coordinates to 1D, get integer index
            if(i == 0)
                index += int(x / _compartmentSize[0]);
            else if(i == 1)
                index += int(x / _compartmentSize[1]) * _grid[0];
            else
                index += int(x / _compartmentSize[2]) * _grid[0] * _grid[1];
            i++;
        }
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
};


#endif /* defined(__Cyto__GController__) */
