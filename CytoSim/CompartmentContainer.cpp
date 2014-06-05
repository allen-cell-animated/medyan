//
//  CompartmentContainer.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 9/5/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "CompartmentContainer.h"

namespace chem {

    template <size_t NDIM>
    void CompartmentsSimpleGrid<NDIM>::generateConnections()
    {
        //Three dimensional
        if (NDIM == 3) {
        
            for(size_t i=0U; i<_grid[0]; ++i) {
        
                for(size_t j=0U; j<_grid[1]; ++j) {
                    
                    for(size_t k=0U; k<_grid[2]; ++k)
                    {
                        Compartment *target = this->getCompartment(i,j,k);
                        
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
        else if (NDIM == 2) {
            
            for(size_t i=0U; i<_grid[0]; ++i) {
                
                for(size_t j=0U; j<_grid[1]; ++j) {
                    
                    Compartment *target = this->getCompartment(i,j);
                    
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
    
    template void CompartmentsSimpleGrid<3>::generateConnections();
    template void CompartmentsSimpleGrid<2>::generateConnections();
    template void CompartmentsSimpleGrid<1>::generateConnections();
    
} // end of chem