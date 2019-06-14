
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "GController.h"

#include "Parser.h"
#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "BoundaryImpl.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"

using namespace mathfunc;
unsigned int GController::getCompartmentID(const vector<floatingpoint> &coords)
{
    //Check if out of bounds
    unsigned int index = 0;
    unsigned int i = 0;
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

    return index;
}

Compartment* GController::getCompartment(const int index){
    try {
        return _compartmentGrid->getCompartment(index);
    }
    catch (exception& e){
        cout << "Bad compartment access at..." << endl;
        cout << "Compartment index = " << index << endl;
//        cout << "Coords = " << coords[0] << " " << coords[1] << " " << coords[2] << endl;
        throw NaNCoordinateException();
    }
}

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

Compartment* GController::getCompartment(const vector<floatingpoint> &coords)
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

void GController::generateConnections()
{
    for(size_t i=0U; i<_grid[0]; ++i) {

        for(size_t j=0U; j<_grid[1]; ++j) {

            for(size_t k=0U; k<_grid[2]; ++k)
            {
                vector<size_t> indices{i,j,k};
                Compartment *target = getCompartment(indices);
                
                vector<floatingpoint> coordinates =
                   {indices[0] * _compartmentSize[0] + _compartmentSize[0] / 2,
                    indices[1] * _compartmentSize[1] + _compartmentSize[1] / 2,
                    indices[2] * _compartmentSize[2] + _compartmentSize[2] / 2};
                target->setCoordinates(coordinates);
                //Go through all 27 neighbors
                int stencilcount = 0;
                for(int ii: {-1,0,1}){
                    for(int jj: {-1,0,1}){
                        for(int kk: {-1,0,1}){
                            //Consider the target bin itself as a neighbor.
                            stencilcount++;
                            int iprime = i+ii;
                            int jprime = j+jj;
                            int kprime = k+kk;

                            if(iprime<0 or iprime==int(_grid[0]) or jprime<0 or
                               jprime==int(_grid[1]) or kprime<0 or
                               kprime==int(_grid[2]))
                                continue;
                            vector<size_t> currentIndices{size_t(iprime), size_t
                                    (jprime), size_t(kprime)};
                            Compartment *neighbor = getCompartment(currentIndices);
                            //27 enclosing neighbors
                            target->addenclosingNeighbour(neighbor, stencilcount -1);
                            if(ii ==0 && jj == 0 && kk ==0) continue;
                            if(jj==0 && kk>=0 && ii <=0)
                                target->adduniquepermuteNeighbour(neighbor, stencilcount -1);
                            else if(jj==0 && kk==1 && ii == 1)
                                target->adduniquepermuteNeighbour(neighbor, stencilcount -1);
                            else if(jj == 1)
                                target->adduniquepermuteNeighbour(neighbor, stencilcount -1);

                            if(ii != 0 && jprime == j && kprime == k)
                                target->addNeighbour(neighbor, (ii < 0? 0: 1));
                            else if(jj != 0 && iprime == i && kprime == k)
                                target->addNeighbour(neighbor, (jj < 0? 2: 3));
                            else if(kk != 0 && iprime == i && jprime == j)
                                target->addNeighbour(neighbor, (kk < 0? 4: 5));
                        }
                    }
                }
                /*for(int ii: {-1,1})
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
                    target->addNeighbour(neighbor);
                }*/
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
                   _compartmentSize[1] * _grid[1] / 2,
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
    vector<BoundaryMove> move;
    for(auto bm:BTypes.boundaryMove){
        if(bm == "NONE") move.push_back(BoundaryMove::None);
        else if(bm == "LEFT") {

#ifndef CHEMISTRY
            cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::Left);
        }
        else if(bm == "RIGHT") {

#ifndef CHEMISTRY
            cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::Right);
        }
        else if(bm == "FRONT") {

#ifndef CHEMISTRY
            cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::Front);
        }
        else if(bm == "BACK") {

#ifndef CHEMISTRY
            cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::Back);
        }
        else if(bm == "BOTTOM") {

#ifndef CHEMISTRY
            cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::Bottom);
        }
        else if(bm == "TOP") {

#ifndef CHEMISTRY
            cout << "Top moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::Top);
        }
        else if(bm == "ALL") {

#ifndef CHEMISTRY
            cout << "Full moving boundary cannot be executed without "
             << "chemistry enabled. Fix these compilation macros "
             << "and try again." << endl;
        exit(EXIT_FAILURE);
#endif
            move.push_back(BoundaryMove::All);
        }
            //if nothing is specified, don't move boundaries
        else if(bm == "") {
            move.push_back(BoundaryMove::None);
        }
        else {
            cout << "Given boundary movement not yet implemented. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
    }

    if(BTypes.boundaryShape == "CUBIC")
        _boundary = new BoundaryCubic(_subSystem, move);

    else if(BTypes.boundaryShape == "SPHERICAL") {

        if(move.size() > 0) {
            if(move[0] != BoundaryMove::None){

                cout << "Moving boundaries for a spherical shape "
                     << "not yet implemented. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }

        _boundary = new BoundarySpherical(_subSystem,
                    SysParams::Boundaries().diameter, move);
    }

    else if(BTypes.boundaryShape == "CAPSULE") {

        if(move.size() > 0) {
            if(move[0] != BoundaryMove::None){

                cout << "Moving boundaries for a capsule shape "
                     << "not yet implemented. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        _boundary = new BoundaryCapsule(_subSystem,
                    SysParams::Boundaries().diameter, move);
    }

    else if(BTypes.boundaryShape == "CYLINDER") {

        if(move.size() > 0) {
            if(move[0] != BoundaryMove::None){

                cout << "Moving boundaries for a cylinder shape "
                << "not yet implemented. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        _boundary = new BoundaryCylinder(_subSystem,
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
        if(_boundary->within(C)) C->setAsActive();
}

void GController::findCompartments(const vector<floatingpoint>& coords,
                                   Compartment* ccheck, floatingpoint dist,
                                   vector<Compartment*>& compartments) {

    //base case : if c and ccheck are not within range, return
    if(twoPointDistancesquared(coords, ccheck->coordinates()) > (dist * dist) ) return;

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
        vector<floatingpoint> coord =
        {_grid[0] * _compartmentSize[0] * Rand::randfloatingpoint(0,0.999),
         _grid[1] * _compartmentSize[1] * Rand::randfloatingpoint(0,0.999),
        _grid[2] * _compartmentSize[2] * Rand::randfloatingpoint(0,0.999)};
        
        Compartment* c = getCompartment(coord);
        if(c->isActivated()) return c;
    }
}

vector<floatingpoint> GController::getRandomCoordinates(Compartment* c) {
    
    //get coordinates of compartment
    auto coordsCompartment = c->coordinates();
    vector<floatingpoint> coords;
    coords.push_back(coordsCompartment[0] + _compartmentSize[0] * Rand::randfloatingpoint(-1,1) / 2);
    coords.push_back(coordsCompartment[1] + _compartmentSize[1] * Rand::randfloatingpoint(-1,1) / 2);
    coords.push_back(coordsCompartment[2] + _compartmentSize[2] * Rand::randfloatingpoint(-1,1) / 2);
    
    return coords;
}

vector<floatingpoint> GController::getRandomCenterCoordinates(Compartment* c) {
    
    //get coordinates of compartment
    auto coordsCompartment = c->coordinates();
    vector<floatingpoint> coords;
    coords.push_back(coordsCompartment[0] + _compartmentSize[0] * Rand::randfloatingpoint(-1,1) / 2);
    coords.push_back(coordsCompartment[1] + _compartmentSize[1] * Rand::randfloatingpoint(-1,1) / 2);
    coords.push_back(coordsCompartment[2] + _compartmentSize[2] * Rand::randfloatingpoint(-0.4,0.4) / 2);
    
    return coords;
}

vector<floatingpoint> GController::getRandomCoordinates() {

    vector<floatingpoint> coords;
    auto bboundsinit = SysParams::Boundaries().fraccompartmentspan;
    coords.push_back(Rand::randfloatingpoint(bboundsinit[0][0],
                                      bboundsinit[1][0]) * _grid[0] * _compartmentSize[0]);
    coords.push_back(Rand::randfloatingpoint(bboundsinit[0][1],
                                      bboundsinit[1][1]) * _grid[1] * _compartmentSize[1]);
    coords.push_back(Rand::randfloatingpoint(bboundsinit[0][2],
                                      bboundsinit[1][2]) * _grid[2] * _compartmentSize[2]);

    return coords;
}

//Qin
vector<floatingpoint> GController::getRandomCenterCoordinates() {

    vector<floatingpoint> coords;

    coords.push_back(Rand::randfloatingpoint(0,1) * _grid[0] * _compartmentSize[0]);
    coords.push_back(Rand::randfloatingpoint(0,1) * _grid[1] * _compartmentSize[1]);
    coords.push_back(Rand::randfloatingpoint(0.3,0.7) * _grid[2] * _compartmentSize[2]);
    
    return coords;
}

short GController::_nDim = 0;

vector<int>    GController::_grid = {};
vector<floatingpoint> GController::_size = {};
vector<floatingpoint> GController::_compartmentSize = {};
vector<floatingpoint> GController::_centerGrid = {};
floatingpoint         GController::_compartmentVolume = 0;
vector<floatingpoint> GController::_compartmentArea = {};

CompartmentGrid* GController::_compartmentGrid = 0;

