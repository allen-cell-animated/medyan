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
    static void initializeGrid(short nDim, std::vector<int> grid, std::vector<double> compartmentSize);
    
    /// Get compartment from the grid
    /// @param - args, the indices in n-dimensions of the compartment
    template<typename ...Args>
    static Compartment* getCompartment(Args&& ...args)
    {
        size_t index = 0;
        size_t i = _nDim-1;
        for(auto x: {args...})
        {
            index+=x*std::pow(_grid[i],i);
            --i;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
    /// Alternate getter from the grid
    static Compartment* getCompartment(const std::vector<size_t> &indices)
    {
        size_t index = 0;
        size_t i = _nDim-1;
        for(auto x: indices)
        {
            index+=x*std::pow(_grid[i],i);
            --i;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
    /// Get the compartment given a set of coordinates
    static Compartment* getCompartment(const std::vector<float> &coords)
    {
        size_t index = 0;
        size_t i = _nDim-1;
        for(auto x: coords)
        {
            index+=int(x / _compartmentSize[index]) * std::pow(_grid[i],i);
            --i;
        }
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
};


#endif /* defined(__Cyto__GController__) */
