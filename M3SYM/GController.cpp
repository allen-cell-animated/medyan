
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

#include "GController.h"

#include "Boundary.h"
#include "Parser.h"
#include "CompartmentContainer.h"

#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;

Compartment* GController::getCompartment(const vector<size_t> &indices)
{
    size_t index = 0;
    size_t i = 0;
    for(auto x: indices)
    {
        //Flatten the indices to 1D
        if(i == 0) {
            if(x >= _grid[0]) {
                cout << endl;
                cout << "Element coordinate = " << indices[0] << ", "
                                                << indices[1] << ", "
                                                << indices[2] << endl;
                throw OutOfBoundsException();
            }
            index += x;
        }
        else if(i == 1) {
            if(x >= _grid[1]) {
                cout << endl;
                cout << "Element coordinate = " << indices[0] << ", "
                                                << indices[1] << ", "
                                                << indices[2] << endl;
                throw OutOfBoundsException();
            }
            index += x * _grid[0];
        }
        else {
            if(x >= _grid[2]) {
                cout << endl;
                cout << "Element coordinate = " << indices[0] << ", "
                                                << indices[1] << ", "
                                                << indices[2] << endl;
                throw OutOfBoundsException();
            }
            index += x * _grid[0] * _grid[1];
        }
        
        i++;
    }
    return static_cast<Compartment*>(
        CompartmentGrid::instance()->children().at(index).get());
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
                cout << endl;
                cout << "Element coordinate = " << coords[0] << ", "
                                                << coords[1] << ", "
                                                << coords[2] << endl;
                throw OutOfBoundsException();
            }
            index += int(x / _compartmentSize[0]);
        }
        else if(i == 1) {
            if(x < 0 || x >= (_compartmentSize[1] * _grid[1])) {
                cout << endl;
                cout << "Element coordinate = " << coords[0] << ", "
                                                << coords[1] << ", "
                                                << coords[2] << endl;
                throw OutOfBoundsException();
            }
            index += int(x / _compartmentSize[1]) * _grid[0];
        }
        else {
            if(x < 0 || x >= (_compartmentSize[2] * _grid[2])) {
                cout << endl;
                cout << "Element coordinate = " << coords[0] << ", "
                                                << coords[1] << ", "
                                                << coords[2] << endl;
                throw OutOfBoundsException();
            }
            index += int(x / _compartmentSize[2]) * _grid[0] * _grid[1];
        }
        i++;
    }
    return static_cast<Compartment*>(
        CompartmentGrid::instance()->children().at(index).get());
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


void GController::initializeGrid() {
    
    //Initial parameters of system
    _nDim = SysParams::Geometry().nDim;
    
    _compartmentSize = {SysParams::Geometry().compartmentSizeX,
                        SysParams::Geometry().compartmentSizeY,
                        SysParams::Geometry().compartmentSizeZ};
    
    _grid = {SysParams::Geometry().NX,
             SysParams::Geometry().NY,
             SysParams::Geometry().NZ};
    
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
        {_grid[0] * _compartmentSize[0] * randomDouble(0,1),
         _grid[1] * _compartmentSize[1] * randomDouble(0,1),
         _grid[2] * _compartmentSize[2] * randomDouble(0,1)};
        
        Compartment* c = getCompartment(coord);
        if(c->isActivated()) return c;
    }
}

vector<double> GController::getRandomCoordinates(Compartment* c) {
    
    //get coordinates of compartment
    auto coordsCompartment = c->coordinates();
    vector<double> coords;
    
    coords.push_back(coordsCompartment[0] +
    _compartmentSize[0] * randomDouble(-1,1) / 2);
    coords.push_back(coordsCompartment[1] +
    _compartmentSize[1] * randomDouble(-1,1) / 2);
    coords.push_back(coordsCompartment[2] +
    _compartmentSize[2] * randomDouble(-1,1) / 2);
    
    return coords;
}



vector<int> GController::_grid = {};
vector<double> GController::_compartmentSize = {};
short GController::_nDim = 0;

