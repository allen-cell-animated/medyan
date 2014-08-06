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
    int _nDim; ///< number of dimensions in system
    std::vector<int> _grid; ///< grid dimensionality
    std::vector<float> _sideLength; ///< side lengths of a compartment
    
    Boundary* _b; ///< boundary of the system
    
    ///Generate all neighbors lists for each compartment
    void generateConnections();
    
public:
    
    ///initialize the grid based on input parameters
    ///@param nDim - the number of dimensions in this system
    ///@param grid - the number of compartments in each dimension
    ///@param systemSize - the actual size of the system in each dimension
    void initialize(int nDim, std::vector<int> grid, std::vector<float> systemSize) {
        
        _nDim = nDim;
        _grid = grid;
        
        int size = 1;
        int i = 0;
        for(auto x: _grid) {
            _sideLength[i] = systemSize[i] / _grid[i];
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
    Compartment* getCompartment(Args&& ...args)
    {
        size_t index = 0;
        size_t i = _nDim-1;
        for(auto x: {args...})
        {
            index+=x*std::pow(_grid[i],i);
            --i;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children()[index].get());
    }
    
    /// Alternate getter from the grid
    Compartment* getCompartment(const std::vector<size_t> &indices) const
    {
        size_t index = 0;
        size_t i = _nDim-1;
        for(auto x: indices)
        {
            index+=x*std::pow(_grid[i],i);
            --i;
        }
        //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children()[index].get());
    }
    
    /// Get the compartment given a set of coordinates
    Compartment* getCompartment(const std::vector<float> &coords) const
    {
        size_t index = 0;
        size_t i = _nDim-1;
        for(auto x: coords)
        {
            index+=int(x / _sideLength[index]) * std::pow(_grid[i],i);
            --i;
        }
        return static_cast<Compartment*>(CompartmentGrid::Instance(CompartmentGridKey())->children()[index].get());
    }
 
};



#endif /* defined(__Cyto__GController__) */
