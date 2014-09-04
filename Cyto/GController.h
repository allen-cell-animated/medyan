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
#include "MathFunctions.h"

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
    ///Generate all neighbors lists for each compartment
    void generateConnections();
    
public:

    ///initialize the grid based on input parameters
    ///@param nDim - the number of dimensions in this system
    ///@param grid - the number of compartments in each dimension
    ///@param systemSize - the actual size of the system in each dimension
    void initializeGrid() {
        
        //make sure dimensions are same as system and grid dimensions
        //assert(nDim == grid.size() && nDim == systemSize.size());
        
        int size = 1;
        int i = 0;
        for(auto x: GRID) {
            size*=x;
            i++;
        }
        ///Set the instance of this grid with given parameters
        CompartmentGrid::setInstance(CompartmentGridKey(), size);
        
        ///Create connections based on dimensionality
        generateConnections();
    }
    
    /// Get compartment from the grid
    /// @param - args, the indices in n-dimensions of the compartment
    template<typename ...Args>
    static Compartment* getCompartment(Args&& ...args)
    {
        size_t index = 0;
        size_t i = NDIM-1;
        for(auto x: {args...})
        {
            index+=x*std::pow(GRID[i],i);
            --i;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
    /// Alternate getter from the grid
    static Compartment* getCompartment(const std::vector<size_t> &indices)
    {
        size_t index = 0;
        size_t i = NDIM-1;
        for(auto x: indices)
        {
            index+=x*std::pow(GRID[i],i);
            --i;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
    
    /// Get the compartment given a set of coordinates
    static Compartment* getCompartment(const std::vector<float> &coords)
    {
        size_t index = 0;
        size_t i = NDIM-1;
        for(auto x: coords)
        {
            index+=int(x / COMPARTMENT_SIZE[index]) * std::pow(GRID[i],i);
            --i;
        }
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children().at(index).get());
    }
};


#endif /* defined(__Cyto__GController__) */
