//
//  CompartmentContainer.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 9/5/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "CompartmentContainer.h"

CompartmentGrid* CompartmentGrid::_instance = 0;

void CompartmentGrid::setInstance(CompartmentGridKey k, std::initializer_list<size_t> grid)
{
    if(_instance != 0)
        delete _instance;
    _instance = new CompartmentGrid(grid);
}

CompartmentGrid* CompartmentGrid::Instance(CompartmentGridKey k) {
    if(_instance==0)
        _instance = new CompartmentGrid({});
    return _instance;
}


void CompartmentGrid::generateConnections()
{
    std::vector<float> sides = _prototype_compartment.sides();
    
    //Three dimensional
    if (_nDim == 3) {
        
        for(size_t i=0U; i<_grid[0]; ++i) {
    
            for(size_t j=0U; j<_grid[1]; ++j) {
                
                for(size_t k=0U; k<_grid[2]; ++k)
                {
                    Compartment *target = this->getCompartment(i,j,k);
                    std::vector<float> coords{i * sides[0], j * sides[1], k * sides[2]};
                    target->setCoords(coords.begin());
                    
                    for(int ii: {-1,1})
                    {
                        int iprime = i+ii;
                        if(iprime<0 or iprime==int(_grid[0]))
                            continue;
                        Compartment *neighbor = this->getCompartment(size_t(iprime),j,k);
                        target->addNeighbour(neighbor);
                    }
                    for(int jj: {-1,1})
                    {
                        int jprime = j+jj;
                        if(jprime<0 or jprime==int(_grid[1]))
                            continue;
                        Compartment *neighbor = this->getCompartment(i,size_t(jprime),k);
                        target->addNeighbour(neighbor);
                    }
                    for(int kk: {-1,1})
                    {
                        int kprime = k+kk;
                        if(kprime<0 or kprime==int(_grid[2]))
                            continue;
                        Compartment *neighbor = this->getCompartment(i,j,size_t(kprime));
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
                
                Compartment *target = this->getCompartment(i,j);
                std::vector<float> coords{i * sides[0], j * sides[1]};
                target->setCoords(coords.begin());
                
                for(int ii: {-1,1})
                {
                    int iprime = i+ii;
                    if(iprime<0 or iprime==int(_grid[0]))
                        continue;
                    Compartment *neighbor = this->getCompartment(size_t(iprime),j);
                    target->addNeighbour(neighbor);
                }
                for(int jj: {-1,1})
                {
                    int jprime = j+jj;
                    if(jprime<0 or jprime==int(_grid[1]))
                        continue;
                    Compartment *neighbor = this->getCompartment(i,size_t(jprime));
                    target->addNeighbour(neighbor);
                }
            }
        }
    }
    
    //One dimensional
    else {
        for(size_t i=0U; i<_grid[0]; ++i) {
        
            Compartment *target = this->getCompartment(i);
            std::vector<float> coords{i * sides[0]};
            target->setCoords(coords.begin());
            
            for(int ii: {-1,1})
            {
                int iprime = i+ii;
                if(iprime<0 or iprime==int(_grid[0]))
                    continue;
                Compartment *neighbor = this->getCompartment(size_t(iprime));
                target->addNeighbour(neighbor);
            }
        }
    }

}