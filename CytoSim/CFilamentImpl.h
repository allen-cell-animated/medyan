//
//  CFilamentImpl.h
//  CytoSim
//
//  Created by James Komianos on 7/17/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CFilamentImpl__
#define __CytoSim__CFilamentImpl__

#include <iostream>
#include "CFilament.h"


namespace chem {
    
    ///Basic monomer consisting of actin, formin, capping, and a virtual front/back
    class CMonomerBasic : public CMonomer {
        
    private:
        SpeciesFilament* _filament_species[1]; ///array of filament species (just actin for now)
        SpeciesFilament* _end_species[4]; ///<array of "end" species (front, back, formin, capping)
        
    public:
        ///Constructor, initializes species container
        CMonomerBasic(std::vector<SpeciesFilament*> species, Compartment* c) : CMonomer(c)
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
        ~CMonomerBasic()
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
        
        ///Look up a species given a name
        virtual SpeciesFilament* getActin() {return _filament_species[0];}
        
        virtual SpeciesFilament* getFront() {return _end_species[0];}
        
        virtual SpeciesFilament* getBack() {return _end_species[1];}
        
        virtual SpeciesFilament* getFormin() {return _end_species[2];}
        
        virtual SpeciesFilament* getCapping() {return _end_species[3];}
        
        ///Check if this monomer is valid
        virtual bool checkSpecies(int sum)
        {
            return true;
//            int currentSum = 0;
//            for(auto &s : _species)
//                currentSum += s->getN();
//            return currentSum = sum;
        }

        ///Print a species in this filament element
        virtual void print()
        {
            for (auto &s : _filament_species)
                if(s->getN() == 1) std::cout << s->getName().at(0);
            for (auto &s : _end_species)
                if(s->getN() == 1) std::cout << s->getName().at(0);
        }
        
        ///move all species (and associated reactions) in this element to another compartment
        virtual void moveToCompartment(Compartment *c) { 
        }
        
        ///Find active filament species
        ///@note return null if none
        virtual Species* getActiveFilamentSpecies() {
            for (auto &s : _filament_species)
                if(s->getN() == 1) return s;
            return nullptr;
        }
        
        ///Find active end species
        ///@note return null if none
        virtual Species* getActiveEndSpecies() {
            for (auto &s : _end_species)
                if(s->getN() == 1) return s;
            return nullptr;
        }
    };
        
    ///Basic bound consisting of myosin, myosin-actin, and a virtual empty species
    class CBoundBasic : public CBound {
        
    private:
        SpeciesBound* _species[3]; ///<array of contained species
        
    public:
        CBoundBasic(std::vector<SpeciesBound*> species, Compartment* c) : CBound(c)
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
        ~CBoundBasic()
        {
            for (auto &s : _species)
            {
                _compartment->removeInternalReactions(s);
                _compartment->removeSpecies(s);
            }
        }
        
        ///Look up a species given a name
        virtual SpeciesBound* getEmpty() {return _species[0];}
        
        virtual SpeciesBound* getMyosin() {return _species[1];}
        
        virtual SpeciesBound* getMyosinActin() {return _species[2];}
        
        ///Check if this monomer is valid
        virtual bool checkSpecies(int sum)
        {
            int currentSum = 0;
            for(auto &s : _species)
                currentSum += s->getN();
            return currentSum = sum;
        }
        
        
        ///Print a species in this filament element
        virtual void print()
        {
            for (auto &s : _species)
                if(s->getN() == 1) std::cout << s->getName().at(0);
        }

        virtual void moveToCompartment(Compartment* c) {}
    };
    
}





#endif /* defined(__CytoSim__CFilamentImpl__) */
