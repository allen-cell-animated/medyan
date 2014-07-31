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

    /// CompartmentGrid class is a simple n-dimensional grid of CompartmentSpatial objects
    
    /*!
     *  The CompartmentGrid class is a grid of CompartmentSpatial objects, each of them seperately
     *  holding internal and diffusion reactions, species information, as well as spatial information.
     *  This class is n-dimensional, and the dimension is specified at runtime.
     *
     *  All compartments within CompartmentGrid are indexed, and this class is responsible for
     *  assignment of compartment neighbors upon initialization.
     *
     *  The _prototype_compartment is used such that all reactions and species can be added to 
     *  the prototype, and this configuration can be copied to all compartments in the grid.
     *  An example of such is below:
     *
     *  @code
     *  CompartmentGrid g{50,50,50};
     *  Compartment &Cproto = g.getProtoCompartment();
     *  Species *M1 = Cproto.addSpecies("Myosin",1U);
     *  Cproto.setDiffusionRate(M1,2000);
     *  Species *M2 = Cproto.addSpecies("Fascin",6U);
     *  Cproto.setDiffusionRate(M2,2000);
     *
     *  vector<float> sides{100.0,100.0,100.0};
     *  Cproto.setSides(sides.begin());
     *
     *  Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
     *  Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
     *  g.initialize();
     *  @endcode
     *
     *  @note the protocompartment's side length must be specified so that all compartments will
     *  be initialized to the same side length. The protocomparment must be set up BEFORE initalize()
     *  is called.
     */
    class CompartmentGrid : public Composite {
    private:
        Compartment _prototype_compartment; ///< prototype compartment, to be configured before initialization
        std::vector<size_t> _grid; ///< size of this grid in all dimensions, in units of compartments
        bool _is_initialized; /// grid is initalized or not
        
        const short _nDim; ///Dimensions of grid
    public:
        
        ///Default constructor
        CompartmentGrid(std::initializer_list<size_t> grid) : _nDim(grid.size()), _prototype_compartment(grid.size()),
                                                              _grid(grid), _is_initialized(false)
        {
            assert(0 < _nDim && _nDim <= 3);
        }
        
        /// Get name of this compartment grid
        virtual std::string getFullName() const {return std::string("CompartmentGrid<")+std::to_string(_nDim)+">";};
        
        /// Initialize the compartment, which copies all species and reactions of the protocompartment into the
        /// compartments in the grid. Also generates neighboring connections for all compartments as well as
        /// initializes spatial coordinates of all compartments.
        /// @note - _prototype_compartment must be configured before this is called
        void initialize()
        {
            if(_is_initialized)
                throw std::runtime_error("CompartmentGrid::initialize() should be called only once");
            size_t length = 1;
            for(auto x: _grid)
            {
                length*=x;
            }
            
            for(size_t i=0; i<length; ++i)
            {
                addChild(std::unique_ptr<Component>(new Compartment(_nDim)));
            }
            
            _is_initialized=true;

            generateConnections();
            
            for(auto &c : children())
            {
                Compartment *C = static_cast<Compartment*>(c.get());
                *C = _prototype_compartment;
            }
        }
        
        ///Activate all compartments
        void activateAll()
        {
            for (auto &c : children())
                static_cast<Compartment*>(c.get())->activate();
            
        }
        
        ///Generate all diffusion reactions for this compartment grid
        ///@note - only compartments that are activated will generate diffusion reactions with one another
        void generateDiffusionReactions()
        {
            for(auto &c : children())
            {
                Compartment *C = static_cast<Compartment*>(c.get());
                C->generateAllDiffusionReactions();
            }
        }
        
        /// Get compartment from the grid
        /// @param - args, the indices in n-dimensions of the compartment
        template<typename ...Args>
        Compartment* getCompartment(Args&& ...args)
        {
            if(not _is_initialized)
                throw std::runtime_error("CompartmentGrid::getCompartment(): initialize() needs to be called first");

            size_t index = 0;
            size_t i = _nDim-1;
            for(auto x: {args...})
            {
                index+=x*std::pow(_grid[i],i);
                --i;
            }
            //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
            return static_cast<Compartment*>(children()[index].get());
        }
        
        /// Alternate getter from the grid
        Compartment* getCompartment(const std::vector<size_t> &indices) const
        {
            if(not _is_initialized)
                throw std::runtime_error("CompartmentGrid::getCompartment(): initialize() needs to be called first");
            
            size_t index = 0;
            size_t i = _nDim-1;
            for(auto x: indices)
            {
                index+=x*std::pow(_grid[i],i);
                --i;
            }
            //            std::cout << "CompartmentGrid::getCompartment(): index=" << index << std::endl;
            return static_cast<Compartment*>(children()[index].get());
        }
        
        /// Get the compartment given a set of coordinates
        Compartment* getCompartment(const std::vector<float> &coords) const
        {
            if(not _is_initialized)
                throw std::runtime_error("CompartmentGrid::getCompartment(): initialize() needs to be called first");
            
            size_t index = 0;
            size_t i = _nDim-1;
            for(auto x: coords)
            {
                index+=(x / getProtoCompartment().sides()[index]) * std::pow(_grid[i],i);
                --i;
            }
            return static_cast<Compartment*>(children()[index].get());
        }
        
        
        /// Get the protocompartment from this grid, in order to configure and then initialize
        Compartment& getProtoCompartment() {return _prototype_compartment;}
        const Compartment& getProtoCompartment() const {return _prototype_compartment;}
        
        /// Generate neighbors and spatial coordinates for all compartments in grid
        void generateConnections();
        
        /// Add reactions to all compartments in the grid
        /// @param - chem, a ChemSim object that controls the reaction algorithm
        virtual void addChemSimReactions()
        {
            for(auto &c : children())
            {
                Compartment*C = static_cast<Compartment*>(c.get());
                C->addChemSimReactions();
            }
        }
        
        /// Print properties of this grid
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
