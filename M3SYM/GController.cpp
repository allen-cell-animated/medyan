
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

#include "GController.h"

#include "Boundary.h"
#include "Parser.h"
#include "CompartmentContainer.h"

#include "MathFunctions.h"
#include "SystemParameters.h"

short GController::_nDim = 0;
vector<int> GController::_grid = {};
vector<double> GController::_compartmentSize = {};


using namespace mathfunc;

Compartment* GController::getCompartment(const vector<size_t> &indices)
{
    size_t index = 0;
    size_t i = 0;
    for(auto x: indices)
    {
        //Flatten the indices to 1D
        if(i == 0) {
            if(x >= _grid[0]) { throw OutOfBoundsException();}
            index += x;
        }
        else if(i == 1) {
            if(x >= _grid[1]) { throw OutOfBoundsException();}
            index += x * _grid[0];
        }
        else {
            if(x >= _grid[2]) { throw OutOfBoundsException();}
            index += x * _grid[0] * _grid[1];
        }
        
        i++;
    }
    return static_cast<Compartment*>(CompartmentGrid::instance()->children().at(index).get());
}

Compartment* GController::getCompartment(const vector<double> &coords)
{
    //Check if out of bounds
    size_t index = 0;
    size_t i = 0;
    for(auto x: coords)
    {
        //Flatten the coordinates to 1D, get integer index
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
    return static_cast<Compartment*>(CompartmentGrid::instance()->children().at(index).get());
}

void GController::generateConnections()
{

    //Three dimensional
    if (_nDim == 3) {
        
        for(size_t i=0U; i<_grid[0]; ++i) {
            
            for(size_t j=0U; j<_grid[1]; ++j) {
                
                for(size_t k=0U; k<_grid[2]; ++k)
                {
                    vector<size_t> indices{i,j,k};
                    Compartment *target = getCompartment(indices);
                    
                    vector<double> coordinates =
                       {indices[0] * _compartmentSize[0] + _compartmentSize[0] / 2,
                        indices[1] * _compartmentSize[1] + _compartmentSize[1] / 2,
                        indices[2] * _compartmentSize[2] + _compartmentSize[2] / 2};
                    target->setCoordinates(coordinates);
                    
                    for(int ii: {-1,1})
                    {
                        int iprime = i+ii;
                        if(iprime<0 or iprime==int(_grid[0]))
                            continue;
                        vector<size_t> currentIndices{size_t(iprime), j, k};
                        Compartment *neighbor = getCompartment(currentIndices);
                        target->addNeighbour(neighbor);
                    }
                    for(int jj: {-1,1})
                    {
                        int jprime = j+jj;
                        if(jprime<0 or jprime==int(_grid[1]))
                            continue;
                        vector<size_t> currentIndices{i, size_t(jprime), k};
                        Compartment *neighbor = getCompartment(currentIndices);
                        target->addNeighbour(neighbor);
                    }
                    for(int kk: {-1,1})
                    {
                        int kprime = k+kk;
                        if(kprime<0 or kprime==int(_grid[2]))
                            continue;
                        vector<size_t> currentIndices{i, j, size_t(kprime)};
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
                
                vector<size_t> indices{i,j};
                Compartment *target = getCompartment(indices);
                
                vector<double> coordinates =
                   {indices[0] * _compartmentSize[0] + _compartmentSize[0] / 2,
                    indices[1] * _compartmentSize[1] + _compartmentSize[1] / 2};
                target->setCoordinates(coordinates);
                
                
                for(int ii: {-1,1})
                {
                    int iprime = i+ii;
                    if(iprime<0 or iprime==int(_grid[0]))
                        continue;
                    vector<size_t> currentIndices{size_t(i), j};
                    Compartment *neighbor = getCompartment(currentIndices);
                    target->addNeighbour(neighbor);
                }
                
                for(int jj: {-1,1})
                {
                    int jprime = j+jj;
                    if(jprime<0 or jprime==int(_grid[1]))
                        continue;
                    vector<size_t> currentIndices{i, size_t(j)};
                    Compartment *neighbor = getCompartment(currentIndices);
                    target->addNeighbour(neighbor);
                }
            }
        }
    }
    
    //One dimensional
    else {
        for(size_t i=0U; i<_grid[0]; ++i) {
            
            vector<size_t> indices{i};
            Compartment *target = getCompartment(indices);
            
            vector<double> coordinates =
                {indices[0] * _compartmentSize[0] + _compartmentSize[0] / 2};
            target->setCoordinates(coordinates);
            
            for(int ii: {-1,1})
            {
                int iprime = i+ii;
                if(iprime<0 or iprime==int(_grid[0]))
                    continue;
                vector<size_t> currentIndices{size_t(i)};
                Compartment *neighbor = getCompartment(currentIndices);
                target->addNeighbour(neighbor);
            }
        }
    }
}


void GController::initializeGrid() {
    
    //Initial parameters of system
    _nDim = SystemParameters::Geometry().nDim;
    
    _grid = {SystemParameters::Geometry().NX,
                             SystemParameters::Geometry().NY,
                             SystemParameters::Geometry().NZ};
    
    _compartmentSize = {SystemParameters::Geometry().compartmentSizeX,
                        SystemParameters::Geometry().compartmentSizeY,
                        SystemParameters::Geometry().compartmentSizeZ};
    
    //Check that grid and compartmentSize match nDim
    if((_nDim == 1 && _grid[0] != 0 && _grid[1] == 0 && _grid[2]==0 &&
        _compartmentSize[0] != 0 && _compartmentSize[1] == 0 && _compartmentSize[2] == 0) ||
       (_nDim == 2 && _grid[0] != 0 && _grid[1] != 0 && _grid[2]==0 &&
        _compartmentSize[0] != 0 && _compartmentSize[1] != 0 && _compartmentSize[2] == 0) ||
       (_nDim == 3 && _grid[0] != 0 && _grid[1] != 0 && _grid[2]!=0 &&
        _compartmentSize[0] != 0 && _compartmentSize[1] != 0 && _compartmentSize[2] != 0)){
    }
    else {
        cout << "Grid parameters are invalid. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    
    int size = 1;
    for(auto x: _grid) {
        if(x != 0) size*=x;
    }
    
    //Set the instance of this grid with given parameters
    CompartmentGrid::setInstance(size);
    
    //Create connections based on dimensionality
    generateConnections();

}

void GController::activateCompartments(Boundary* boundary) {
    
    //initialize all compartments equivalent to cproto
    for(auto &c : CompartmentGrid::instance()->children()) {
        Compartment *C = static_cast<Compartment*>(c.get());
        if(boundary->within(C->coordinates())) C->activate();
    }
}

void GController::findCompartments(const vector<double>& coords, Compartment* ccheck,
                                   double dist, vector<Compartment*>& compartments) {
    
    //base case : if c and ccheck are not within range, return
    if(twoPointDistance(coords, ccheck->coordinates()) > dist ) return;
    
    //recursive case, c and ccheck are in range. call for all neighbors
    else {
        //if not already in list, add it
        auto it = find(compartments.begin(), compartments.end(), ccheck);
        
        if( it == compartments.end()) {
            //add the compartment
            compartments.push_back(ccheck);
            
            //recursively call for all neighbors
            for(auto &n : ccheck->getNeighbours())
                findCompartments(coords, n, dist, compartments);
        }
    }
}

