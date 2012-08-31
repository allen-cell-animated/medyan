//
//  SpeciesContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 8/30/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__SpeciesContainer__
#define __CytoSim__SpeciesContainer__

#include <iostream>
#include "Species.h"

namespace  chem {

    class SpeciesContainerIFace {
    public:
        virtual void removeSpecies(const std::string &name) = 0;
        virtual Species& findSpecies(const std::string &name) = 0;
        virtual size_t findSpeciesIndex(const std::string &name) = 0;
        virtual Species& findSpecies (size_t index) = 0;
        virtual void printSpecies() {}
    };
    
    class SpeciesPtrContainerIFace {
    public:
        virtual Species* addSpeciesUnique(std::unique_ptr<Species> &&species) = 0;
//        virtual Species* addSpecies(const std::string &name, species_copy_t copy) = 0;
        virtual void removeSpecies(const std::string &name) = 0;
        virtual void removeSpecies(Species* species) = 0;
        virtual Species* findSpecies(const std::string &name) = 0;
        virtual Species* findSpecies (size_t index) = 0;
        virtual void printSpecies() {}
    };
    
    
    

    class SpeciesPtrContainerVector : public  SpeciesPtrContainerIFace {
    protected:
        std::vector<std::unique_ptr<Species>> _species;
    public:
        virtual Species* addSpeciesUnique (std::unique_ptr<Species> &&species) {
            _species.push_back(std::move(species));
            return _species.back().get();
        }
        
        template<typename T, typename ...Args>
        Species* addSpecies( Args&& ...args )
        {
            _species.push_back(std::unique_ptr<T>( new T( std::forward<Args>(args)...) ));
            //        _species.emplace_back(make_unique(Args...));
            return _species.back().get();
        }
        
        virtual void removeSpecies(const std::string &name) {
            auto child_iter = std::find_if(_species.begin(),_species.end(),
                                           [&name](const std::unique_ptr<Species> &element)
                                           {
                                               return element->getName()==name ? true : false;
                                           });
            if(child_iter!=_species.end())
                _species.erase(child_iter);
        }
        
        virtual void removeSpecies (Species* species) {
            auto child_iter = std::find_if(_species.begin(),_species.end(),
                                           [&species](const std::unique_ptr<Species> &element)
                                           {
                                               return element.get()==species ? true : false;
                                           });
            if(child_iter!=_species.end())
                _species.erase(child_iter);
        }
        
        virtual Species* findSpecies(const std::string &name) {
            auto child_iter = std::find_if(_species.begin(),_species.end(),
                                           [&name](const std::unique_ptr<Species> &element)
                                           {
                                               return element->getName()==name ? true : false;
                                           });
            if(child_iter!=_species.end())
                return child_iter->get();
            else
                return nullptr;
        }
        
        virtual Species* findSpecies (size_t index) {
            return _species[index].get();
        }
        
        virtual void printSpecies() {
            for(auto &s : _species)
                std::cout << (*s.get()) << std::endl;
        }

    };
    
    template <class SpeciesSpecific>
    class SpeciesContainerVector : public  SpeciesContainerIFace {
    protected:
        std::vector<SpeciesSpecific> _species;
    public:
        template<typename ...Args>
        size_t addSpecies(Args&& ...args){
//            std::cout << "SpeciesContainerVector::addSpecies()..." << std::endl;
            _species.push_back({std::forward<Args>(args)...});
            return _species.size()-1;
        }
        
        virtual void removeSpecies(const std::string &name) {
            auto child_iter = std::find_if(_species.begin(),_species.end(),
                                           [&name](const Species &element)
                                           {
                                               return element.getName()==name ? true : false;
                                           });
            if(child_iter!=_species.end())
                _species.erase(child_iter);
        }
        
        virtual Species& findSpecies(const std::string &name) {
            auto child_iter = std::find_if(_species.begin(),_species.end(),
                                           [&name](const Species &element)
                                           {
                                               return element.getName()==name ? true : false;
                                           });
            if(child_iter!=_species.end())
                return (*child_iter);
            else
                throw std::out_of_range("Species::findSpecies(): The name was not found");
        }
        
        virtual size_t findSpeciesIndex(const std::string &name) {
            size_t index = 0;
            for(auto &s : _species){
                if(s.getName()==name)
                    return index;
                else
                    ++index;
            }
            throw std::out_of_range("Species::findSpecies(): The name was not found");
        }

        
        virtual Species& findSpecies (size_t index) {
            return _species[index];
        }
        
        virtual void printSpecies() {
            for(auto &s : _species)
                std::cout << s << std::endl;
        }

    };
    
} // end of chem

#endif /* defined(__CytoSim__SpeciesContainer__) */
