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
#include "ChemSim.h"

namespace chem {

    template <size_t NDIM>
    class CompartmentsSimpleGrid : public Composite {
    private:
        CompartmentNDim<NDIM> _prototype_compartment;
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
            
            //TODO
            //Initialize coords and sizes here!!!!
            
            
            for(size_t i=0; i<length; ++i)
            {
                addChild(std::unique_ptr<Component>(new Compartment()));
            }
            
            _is_initialized=true;

            generateConnections();
            
            for(auto &c : children())
            {
                Compartment *C = static_cast<Compartment*>(c.get());
                *C = _prototype_compartment;
            }
            
            for(auto &c : children())
            {
                Compartment *C = static_cast<Compartment*>(c.get());
                C->generateDiffusionReactions();
            }
            
//           std::cout << "CompartmentsSimpleGrid(CTOR): length=" << length << std::endl;
            
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
                index+=x*std::pow(_grid[i],i);
                --i;
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
                index+=x*std::pow(_grid[i],i);
                --i;
            }
            //            std::cout << "CompartmentsSimpleGrid::getCompartment(): index=" << index << std::endl;
            return static_cast<Compartment*>(children()[index].get());
        }
        
        Compartment& getProtoCompartment() {return _prototype_compartment;}
        const Compartment& getProtoCompartment() const {return _prototype_compartment;}
        
        void generateConnections();
        
        virtual void addChemSimReactions(ChemSim &chem)
        {
            for(auto &c : children())
            {
                Compartment *C = static_cast<Compartment*>(c.get());
                C->addChemSimReactions(chem);
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
