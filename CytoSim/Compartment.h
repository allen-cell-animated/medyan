//
//  Compartment.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Compartment_h
#define CytoSim_Experimenting_Compartment_h

#include <vector>
#include "Reaction.h"
#include "Species.h"
#include "SpeciesContainer.h"
#include "Composite.h"

namespace chem {
    
    class Compartment : public Composite, public SpeciesContainerVector<SpeciesDiffusing> {
    protected:
        std::vector<Compartment*> _neighbours;
    public:
        Compartment() :  Composite(), SpeciesContainerVector<SpeciesDiffusing>() {}
        virtual ~Compartment() noexcept {}
        
        virtual bool isSpeciesContainer() const {return true;}
        virtual bool isReactionsContainer() const {return true;}
        virtual std::string getFullName() const {return "Compartment";};
        virtual size_t numberOfSpecies() const {return _species.size();}
        virtual size_t numberOfReactions() const {return 0;}
        
        template<typename ...Args>
        size_t addSpecies(Args&& ...args){
//            std::cout << "Compartment::addSpecies()..." << std::endl;
            size_t index = SpeciesContainerVector<SpeciesDiffusing>::addSpecies(std::forward<Args>(args)...);
            _species[index].setParent(this);
            return index;
        }
        
        std::vector<SpeciesDiffusing>& species() {return _species;}
        const std::vector<SpeciesDiffusing>& species() const {return _species;}
        
        void addNeighbour(Compartment *comp) {
            auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
            if(nit==_neighbours.end())
                _neighbours.push_back(comp);
            else
                throw std::runtime_error("Compartment::addNeighbour(): Compartment is already a neighbour");
        }
        
        void removeNeighbour(Compartment *comp) {
            auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
            if(nit!=_neighbours.end())
                _neighbours.erase(nit);
            else
                throw std::out_of_range("Compartment::removeNeighbour(): Compartment is not a neighbour");
        }
        
        size_t numberOfNeighbours() const {return _neighbours.size();}
        
        std::vector<Compartment*>& neighbours() {return _neighbours;}
        
        const std::vector<Compartment*>& neighbours() const {return _neighbours;}
        
    };



}// end of chem
//

#endif
