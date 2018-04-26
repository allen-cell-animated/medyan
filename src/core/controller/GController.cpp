
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "core/controller/GController.h"

#include "Parser.h"
#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "BoundaryImpl.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"

using namespace mathfunc;

Compartment* GController::getCompartment(const vector<size_t> &indices)
{
    size_t index = 0;
    size_t i = 0;
    for(auto x: indices)
    {
        //Flatten the indices to 1D
        if(i == 0) {
            if(x >= _grid[0])
                throw OutOfBoundsException();

            index += x;
        }
        else if(i == 1) {
            if(x >= _grid[1])
                throw OutOfBoundsException();
        
            index += x * _grid[0];
        }
        else {
            if(x >= _grid[2])
                throw OutOfBoundsException();
            
            index += x * _grid[0] * _grid[1];
        }
        
        i++;
    }
    try {
        return _compartmentGrid->getCompartment(index);
    }
    catch (exception& e){
        cout << "Bad compartment access at..." << endl;
        cout << "Compartment index = " << index << endl;
        cout << "Indices = " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
        throw NaNCoordinateException();
    }
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
            if(x < 0 || x >= (_compartmentSize[0] * _grid[0]))
                throw OutOfBoundsException();
            
            index += int(x / _compartmentSize[0]);
        }
        else if(i == 1) {
            if(x < 0 || x >= (_compartmentSize[1] * _grid[1]))
                throw OutOfBoundsException();
            
            index += int(x / _compartmentSize[1]) * _grid[0];
        }
        else {
            if(x < 0 || x >= (_compartmentSize[2] * _grid[2]))
                throw OutOfBoundsException();
        
            index += int(x / _compartmentSize[2]) * _grid[0] * _grid[1];
        }
        i++;
    }
    
    try {
        return _compartmentGrid->getCompartment(index);
    }
    catch (exception& e){
        cout << "Bad compartment access at..." << endl;
        cout << "Compartment index = " << index << endl;
        cout << "Coords = " << coords[0] << " " << coords[1] << " " << coords[2] << endl;
        throw NaNCoordinateException();
    }
}

vector<size_t> GController::getCompartmentIndices(const vector<double>& coords) {
    /**************************************************************************
    The function takes the coordinates of a point and return the compartment
    indices. All related variables must be with dimension _nDim.

    This function checks boundary.
    **************************************************************************/
    // Assert coords.size() == _grid.size() == _compartmentSize.size() == _nDim
    vector<size_t> indices(_nDim);
    for(size_t coordIdx = 0; coordIdx < _nDim; ++coordIdx) {
        if(coords[coordIdx] < 0 || coords[coordIdx] >= (_compartmentSize[coordIdx] * _grid[coordIdx]))
            throw OutOfBoundsException();
            
        indices[coordIdx] = size_t(int(coords[coordIdx] / _compartmentSize[coordIdx]));
    }
    return indices;
}
bool GController::indicesOutOfBound(const vector<size_t>& indices) {
    for(size_t coordIdx = 0; coordIdx < _nDim; ++coordIdx)
        if(indices[coordIdx] >= _grid[coordIdx])
            return true;
    return false;
}

void GController::generateConnections()
{
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
                    target->addNeighbour(neighbor, (ii < 0? 0: 1));
                }
                for(int jj: {-1,1})
                {
                    int jprime = j+jj;
                    if(jprime<0 or jprime==int(_grid[1]))
                        continue;
                    vector<size_t> currentIndices{i, size_t(jprime), k};
                    Compartment *neighbor = getCompartment(currentIndices);
                    target->addNeighbour(neighbor, (jj < 0? 2: 3));
                }
                for(int kk: {-1,1})
                {
                    int kprime = k+kk;
                    if(kprime<0 or kprime==int(_grid[2]))
                        continue;
                    vector<size_t> currentIndices{i, j, size_t(kprime)};
                    Compartment *neighbor = getCompartment(currentIndices);
                    target->addNeighbour(neighbor, (kk < 0? 4: 5));
                }
            }
        }
    }
    
}


CompartmentGrid* GController::initializeGrid() {
    
    //Initial parameters of system
    _nDim = SysParams::Geometry().nDim;
    
    _compartmentSize = {SysParams::Geometry().compartmentSizeX,
                        SysParams::Geometry().compartmentSizeY,
                        SysParams::Geometry().compartmentSizeZ};
    
    _grid = {SysParams::Geometry().NX,
             SysParams::Geometry().NY,
             SysParams::Geometry().NZ};
    
    _size = {_compartmentSize[0] * _grid[0],
             _compartmentSize[1] * _grid[1],
             _compartmentSize[2] * _grid[2]};
    
    _centerGrid = {_compartmentSize[0] * _grid[0] / 2,
                   _compartmentSize[1] * _grid[2] / 2,
                   _compartmentSize[2] * _grid[2] / 2};
    
    _compartmentVolume = _compartmentSize[0] * _compartmentSize[1] * _compartmentSize[2];

    _compartmentArea = {{
        _compartmentSize[1] * _compartmentSize[2],
        _compartmentSize[2] * _compartmentSize[0],
        _compartmentSize[0] * _compartmentSize[1]
    }};
    
    //Check that grid and compartmentSize match nDim
    if((_nDim == 3 &&
        _grid[0] != 0 && _grid[1] != 0 && _grid[2]!=0 &&
        _compartmentSize[0] != 0 &&
        _compartmentSize[1] != 0 &&
        _compartmentSize[2] != 0)){
    }
    else {
        cout << "Grid parameters are invalid. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    int size = 1;
    for(auto x: _grid) {
        if(x != 0) size*=x;
    }
    
    //Set the instance of this grid with given parameters
    _compartmentGrid = new CompartmentGrid(size);
    
    //Create connections based on dimensionality
    generateConnections();
    
    _subSystem->setCompartmentGrid(_compartmentGrid);
    return _compartmentGrid;
}

Boundary* GController::initializeBoundary(BoundaryType& BTypes) {
    
    BoundaryType type;
    BoundaryMove move;
    
    if(BTypes.boundaryMove == "NONE") move = BoundaryMove::None;
    else if(BTypes.boundaryMove == "TOP") {
        
#ifndef CHEMISTRY
        cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
        move = BoundaryMove::Top;
    }
    else if(BTypes.boundaryMove == "ALL") {
        
#ifndef CHEMISTRY
        cout << "Full moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
        
        move = BoundaryMove::All;
    }
    //if nothing is specified, don't move boundaries
    else if(BTypes.boundaryMove == "") {
        move = BoundaryMove::None;
    }
    else {
        cout << "Given boundary movement not yet implemented. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(BTypes.boundaryShape == "CUBIC")
        _boundary = new BoundaryCubic(_subSystem, move);
    
    else if(BTypes.boundaryShape == "SPHERICAL") {
        
        if(move != BoundaryMove::None) {
            
            cout << "Moving boundaries for a spherical shape "
                 << "not yet implemented. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        _boundary = new BoundarySpherical(_subSystem,
                    SysParams::Boundaries().diameter, move);
    }
    
    else if(BTypes.boundaryShape == "CAPSULE") {
        
        if(move != BoundaryMove::None) {
            
            cout << "Moving boundaries for a capsule shape "
                 << "not yet implemented. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        _boundary = new BoundaryCapsule(_subSystem,
                    SysParams::Boundaries().diameter, move);
    }
    else{
        cout << endl << "Given boundary shape not yet implemented. Exiting." <<endl;
        exit(EXIT_FAILURE);
    }
    
    _subSystem->addBoundary(_boundary);
    return _boundary;
}

void GController::setActiveCompartments() {

    //initialize all compartments equivalent to cproto
    for(auto C : _compartmentGrid->getCompartments())
        if(_boundary->within(C))  C->setAsActive();
}

void GController::findCompartments(const vector<double>& coords,
                                   Compartment* ccheck, double dist,
                                   vector<Compartment*>& compartments) {
    
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

Compartment* GController::getRandomCompartment() {
    
    //return a compartment that is activated
    while(true) {
        
        //create a random coordinate
        vector<double> coord =
        {_grid[0] * _compartmentSize[0] * Rand::randDouble(0,0.999),
         _grid[1] * _compartmentSize[1] * Rand::randDouble(0,0.999),
         _grid[2] * _compartmentSize[2] * Rand::randDouble(0,0.999)};
        
        Compartment* c = getCompartment(coord);
        if(c->isActivated()) return c;
    }
}

vector<double> GController::getRandomCoordinates(Compartment* c) {
    
    //get coordinates of compartment
    auto coordsCompartment = c->coordinates();
    vector<double> coords;
    
    coords.push_back(coordsCompartment[0] + _compartmentSize[0] * Rand::randDouble(-1,1) / 2);
    coords.push_back(coordsCompartment[1] + _compartmentSize[1] * Rand::randDouble(-1,1) / 2);
    coords.push_back(coordsCompartment[2] + _compartmentSize[2] * Rand::randDouble(-1,1) / 2);
    
    return coords;
}

vector<double> GController::getRandomCoordinates() {
    
    vector<double> coords;
    
    coords.push_back(Rand::randDouble(0,1) * _grid[0] * _compartmentSize[0]);
    coords.push_back(Rand::randDouble(0,1) * _grid[1] * _compartmentSize[1]);
    coords.push_back(Rand::randDouble(0,1) * _grid[2] * _compartmentSize[2]);
    
    return coords;
}

short GController::_nDim = 0;

vector<int>    GController::_grid = {};
vector<double> GController::_size = {};
vector<double> GController::_compartmentSize = {};
vector<double> GController::_centerGrid = {};
double         GController::_compartmentVolume = 0;
vector<double> GController::_compartmentArea = {};

CompartmentGrid* GController::_compartmentGrid = 0;

