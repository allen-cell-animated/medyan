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
    void CompartmentsSimpleGrid<NDIM>::generateConnections3D()
    {
        assert(NDIM==3);
        for(size_t i=0U; i<_grid[0]; ++i)
            for(size_t j=0U; j<_grid[1]; ++j)
                for(size_t k=0U; k<_grid[2]; ++k)
                {
                    Compartment *target = this->getCompartment(i,j,k);
                    //                        std::cout << "CompartmentsSimpleGrid::generateConnections3D(): " << i << " " << j << " " << k << std::endl;
                    for(int ii: {-1,1})
                    {
                        int iprime = i+ii;
                        if(iprime<0 or iprime==int(_grid[0]))
                            continue;
                        Compartment *neighbor = this->getCompartment(size_t(iprime),j,k);
                        target->addNeighbour(neighbor);
                        //                            std::cout << "Added neighbor, " << iprime << " " << j << " " << k << std::endl;
                    }
                    for(int jj: {-1,1})
                    {
                        int jprime = j+jj;
                        if(jprime<0 or jprime==int(_grid[1]))
                            continue;
                        Compartment *neighbor = this->getCompartment(i,size_t(jprime),k);
                        target->addNeighbour(neighbor);
                        //                            std::cout << "Added neighbor, " << i << " " << jprime << " " << k << std::endl;
                    }
                    for(int kk: {-1,1})
                    {
                        int kprime = k+kk;
                        if(kprime<0 or kprime==int(_grid[2]))
                            continue;
                        Compartment *neighbor = this->getCompartment(i,j,size_t(kprime));
                        target->addNeighbour(neighbor);
                        //                            std::cout << "Added neighbor, " << i << " " << j << " " << kprime << std::endl;
                    }
                }
        
    }

    template void CompartmentsSimpleGrid<3>::generateConnections3D();
    
} // end of chem