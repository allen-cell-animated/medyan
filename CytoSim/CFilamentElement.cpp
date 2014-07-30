//
//  CFilamentElementImpl.cpp
//  CytoSim
//
//  Created by James Komianos on 7/17/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CFilamentElementImpl.h"

namespace chem {
    
    ///Constructor, initializes species container
    CMonomerBasic::CMonomerBasic(std::vector<SpeciesFilament*> species, Compartment* c) : CMonomer(c)
    {
        ///Initialize member array of species
        for(auto &s : species) {
            std::string name = s->getName();
            
            ///Set up array
            if(name == "Actin")
                _filament_species[0] = s;
            
            else if(name == "Front")
                _end_species[0] = s;
            
            else if(name == "Back")
                _end_species[1] = s;
            
            else if(name == "X-Formin")
                _end_species[2] = s;
            
            else if(name == "Capping")
                _end_species[3] = s;
            
            else {}
        }
    }
    ///Destructor, removes all species and associated reactions from compartment
    CMonomerBasic::~CMonomerBasic()
    {
        for (auto &s : _filament_species)
        {
            _compartment->removeInternalReactions(s);
            _compartment->removeSpecies(s);
        }
        
        for(auto &s : _end_species)
        {
            _compartment->removeInternalReactions(s);
            _compartment->removeSpecies(s);
        }
    }
    
    ///Look up species by name
    Species* CMonomerBasic::getSpeciesByName(std::string& name)
    {
        if(name == "Actin")
            return _filament_species[0];
        
        else if(name == "Front")
            return _end_species[0];
        
        else if(name == "Back")
            return _end_species[1];
        
        else if(name == "X-Formin")
            return _end_species[2];
        
        else if(name == "Capping")
            return _end_species[3];
        
        else {return nullptr;}
    }
    
    ///Find active filament species
    ///@note return null if none
    Species* CMonomerBasic::getActiveFilamentSpecies() {
        for (auto &s : _filament_species)
            if(s->getN() == 1) return s;
        return nullptr;
    }
    
    ///Find active end species
    ///@note return null if none
    Species* CMonomerBasic::getActiveEndSpecies() {
        for (auto &s : _end_species)
            if(s->getN() == 1) return s;
        return nullptr;
    }
    
    CBoundBasic::CBoundBasic(std::vector<SpeciesBound*> species, Compartment* c) : CBound(c)
    {
        ///Initialize member array of species
        for(auto &s : species) {
            std::string name = s->getName();
            
            ///Set up array
            if(name == "Empty")
                _species[0] = s;
            
            else if(name == "Myosin")
                _species[1] = s;
            
            else if(name == "A-MyosinActin")
                _species[2] = s;
            
            else {}
            
        }
    }
    
    ///Destructor, removes all species and associated reactions from compartment
    CBoundBasic::~CBoundBasic()
    {
        for (auto &s : _species)
        {
            _compartment->removeInternalReactions(s);
            _compartment->removeSpecies(s);
        }
    }
    
    ///Check if this monomer is valid
    bool CBoundBasic::checkSpecies(int sum)
    {
        int currentSum = 0;
        for(auto &s : _species)
            currentSum += s->getN();
        return currentSum = sum;
    }
    
    
    ///Look up species by name
    Species* CBoundBasic::getSpeciesByName(std::string& name)
    {
        if(name == "Empty")
            return _species[0];
        
        else if(name == "Myosin")
            return _species[1];
        
        else if(name == "A-MyosinActin")
            return _species[2];
        
        else {return nullptr;}
    }
    
    ///Print a species in this filament element
    void CBoundBasic::print()
    {
        for (auto &s : _species)
            if(s->getN() == 1) std::cout << s->getName().at(0);
    }
    
    
    
    
}