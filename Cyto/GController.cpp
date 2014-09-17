//
//  GController.cpp
//  Cyto
//
//  Created by James Komianos on 8/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "GController.h"
#include "Parser.h"

void GController::initializeGrid() {
    
    ///Initial parameters of system
    short nDim = SystemParameters::Geometry().nDim;
    std::vector<int> grid = {SystemParameters::Geometry().NX,
        SystemParameters::Geometry().NY,
        SystemParameters::Geometry().NZ};
    std::vector<double> compartmentSize = {SystemParameters::Geometry().compartmentSizeX,
        SystemParameters::Geometry().compartmentSizeY,
        SystemParameters::Geometry().compartmentSizeZ};
    
    ///Check that grid and compartmentSize match nDim
    if((nDim == 1 && grid[0] != 0 && grid[1] == 0 && grid[2]==0 && compartmentSize[0] != 0 && compartmentSize[1] == 0 && compartmentSize[2] == 0)
           || (nDim == 2 && grid[0] != 0 && grid[1] != 0 && grid[2]==0 && compartmentSize[0] != 0 && compartmentSize[1] != 0 && compartmentSize[2] == 0)
           || (nDim == 3 && grid[0] != 0 && grid[1] != 0 && grid[2]!=0 && compartmentSize[0] != 0 && compartmentSize[1] != 0 && compartmentSize[2] != 0)){
    }
    else {
        std::cout << "Grid parameters are invalid. Exiting" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    _nDim = nDim;
    
    /// set up
    if(_nDim == 1) {
        _grid = {grid[0]};
        _compartmentSize = {compartmentSize[0]};
    }
    else if(_nDim == 2) {
        _grid = {grid[0], grid[1]};
        _compartmentSize = {compartmentSize[0], compartmentSize[1]};
        
    }
    else {
        _grid = {grid[0], grid[1], grid[2]};
        _compartmentSize = {compartmentSize[0], compartmentSize[1], compartmentSize[2]};
    }
    
    int size = 1;
    int i = 0;
    for(auto x: _grid) {
        size*=x;
        i++;
    }
    
    ///Set the instance of this grid with given parameters
    CompartmentGrid::setInstance(CompartmentGridKey(), size);
    
    ///Create connections based on dimensionality
    generateConnections();
}


void GController::generateConnections()
{
    
    //Three dimensional
    if (_nDim == 3) {
        
        for(size_t i=0U; i<_grid[0]; ++i) {
            
            for(size_t j=0U; j<_grid[1]; ++j) {
                
                for(size_t k=0U; k<_grid[2]; ++k)
                {
                    std::vector<size_t> indices{i,j,k};
                    Compartment *target = getCompartment(indices);
                    
                    for(int ii: {-1,1})
                    {
                        int iprime = i+ii;
                        if(iprime<0 or iprime==int(_grid[0]))
                            continue;
                        std::vector<size_t> currentIndices{size_t(iprime), j, k};
                        Compartment *neighbor = getCompartment(currentIndices);
                        target->addNeighbour(neighbor);
                    }
                    for(int jj: {-1,1})
                    {
                        int jprime = j+jj;
                        if(jprime<0 or jprime==int(_grid[1]))
                            continue;
                        std::vector<size_t> currentIndices{i, size_t(jprime), k};
                        Compartment *neighbor = getCompartment(currentIndices);
                        target->addNeighbour(neighbor);
                    }
                    for(int kk: {-1,1})
                    {
                        int kprime = k+kk;
                        if(kprime<0 or kprime==int(_grid[2]))
                            continue;
                        std::vector<size_t> currentIndices{i, j, size_t(kprime)};
                        Compartment *neighbor = getCompartment(currentIndices);
                        target->addNeighbour(neighbor);
                    }
                }
            }
        }
    }
    
    //Two dimensional
    else if (_nDim == 2) {
        
        for(size_t i=0U; i<_grid[0]; ++i) {
            
            for(size_t j=0U; j<_grid[1]; ++j) {
                
                std::vector<size_t> indices{i,j};
                Compartment *target = getCompartment(indices);
                
                for(int ii: {-1,1})
                {
                    int iprime = i+ii;
                    if(iprime<0 or iprime==int(_grid[0]))
                        continue;
                    std::vector<size_t> currentIndices{size_t(i), j};
                    Compartment *neighbor = getCompartment(currentIndices);
                    target->addNeighbour(neighbor);
                }
                
                for(int jj: {-1,1})
                {
                    int jprime = j+jj;
                    if(jprime<0 or jprime==int(_grid[1]))
                        continue;
                    std::vector<size_t> currentIndices{i, size_t(j)};
                    Compartment *neighbor = getCompartment(currentIndices);
                    target->addNeighbour(neighbor);
                }
            }
        }
    }
    
    //One dimensional
    else {
        for(size_t i=0U; i<_grid[0]; ++i) {
            
            std::vector<size_t> indices{i};
            Compartment *target = getCompartment(indices);
            
            for(int ii: {-1,1})
            {
                int iprime = i+ii;
                if(iprime<0 or iprime==int(_grid[0]))
                    continue;
                std::vector<size_t> currentIndices{size_t(i)};
                Compartment *neighbor = getCompartment(currentIndices);
                target->addNeighbour(neighbor);
            }
        }
    }
}

short GController::_nDim = 0;
std::vector<int> GController::_grid = {};
std::vector<double> GController::_compartmentSize = {};

