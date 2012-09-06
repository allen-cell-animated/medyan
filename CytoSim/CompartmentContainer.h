//
//  CompartmentContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 9/5/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CompartmentContainer__
#define __CytoSim__CompartmentContainer__

#include <iostream>
#include "Compartment.h"

namespace chem {

    template <size_t NDIM>
    class CompartmentsSimpleGrid : public Composite {
    private:
        Compartment _prototype_compartment;
        std::vector<size_t> _grid;
        bool _is_initialized;
    public:
        CompartmentsSimpleGrid(std::initializer_list<size_t> grid) : _grid(grid), _is_initialized(false)
        {
            assert(_grid.size()==NDIM);
        }
        
        virtual std::string getFullName() const {return std::string("CompartmentsSimpleGrid<")+std::to_string(NDIM)+">";};
        
        void initialize()
        {
            if(_is_initialized)
                throw std::runtime_error("CompartmentsSimpleGrid::initialize() should be called only once");
            size_t length = 1;
            for(auto x: _grid)
            {
                length*=x;
            }
            
            for(size_t i=0; i<length; ++i)
            {
                addChild(std::unique_ptr<Component>(new Compartment(_prototype_compartment)));
            }
            
//           std::cout << "CompartmentsSimpleGrid(CTOR): length=" << length << std::endl;
            
            _is_initialized=true;

            if(NDIM==3)
                generateConnections3D();
            
            for(auto &c: children())
            {
                Compartment *C = static_cast<Compartment*>(c.get());
                C->generateDiffusionReactions();
            }
        }
        
        template<typename ...Args>
        Compartment* getCompartment(Args&& ...args)
        {
            if(not _is_initialized)
                throw std::runtime_error("CompartmentsSimpleGrid::getCompartment(): initialize() needs to be called first");

            size_t index = 0;
            size_t i = NDIM-1;
            for(auto x: {args...})
            {
                index+=x*std::pow(_grid[i],i--);
            }
            //            std::cout << "CompartmentsSimpleGrid::getCompartment(): index=" << index << std::endl;
            return static_cast<Compartment*>(children()[index].get());
        }
        
        Compartment* getCompartment(const std::vector<size_t> &indices) const
        {
            if(not _is_initialized)
                throw std::runtime_error("CompartmentsSimpleGrid::getCompartment(): initialize() needs to be called first");
            
            size_t index = 0;
            size_t i = NDIM-1;
            for(auto x: indices)
            {
                index+=x*std::pow(_grid[i],i--);
            }
            //            std::cout << "CompartmentsSimpleGrid::getCompartment(): index=" << index << std::endl;
            return static_cast<Compartment*>(children()[index].get());
        }
        
        Compartment& getProtoCompartment() {return _prototype_compartment;}
        const Compartment& getProtoCompartment() const {return _prototype_compartment;}
        
        void generateConnections3D()
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
                            if(iprime<0 or iprime==_grid[0])
                                continue;
                            Compartment *neighbor = this->getCompartment(size_t(iprime),j,k);
                            target->addNeighbour(neighbor);
//                            std::cout << "Added neighbor, " << iprime << " " << j << " " << k << std::endl;
                        }
                        for(int jj: {-1,1})
                        {
                            int jprime = j+jj;
                            if(jprime<0 or jprime==_grid[1])
                                continue;
                            Compartment *neighbor = this->getCompartment(i,size_t(jprime),k);
                            target->addNeighbour(neighbor);
//                            std::cout << "Added neighbor, " << i << " " << jprime << " " << k << std::endl;
                        }
                        for(int kk: {-1,1})
                        {
                            int kprime = k+kk;
                            if(kprime<0 or kprime==_grid[2])
                                continue;
                            Compartment *neighbor = this->getCompartment(i,j,size_t(kprime));
                            target->addNeighbour(neighbor);
//                            std::cout << "Added neighbor, " << i << " " << j << " " << kprime << std::endl;
                        }
                    }

        }
        
        virtual void printSelf()
        {
            std::cout << getFullName() << std::endl;
            std::cout << "Number of Compartment objects: " << numberOfChildren() << std::endl;
            for(auto &c : children())
                c->printSelf();
        }

    };


}// end of chem
//



#endif /* defined(__CytoSim__CompartmentContainer__) */
