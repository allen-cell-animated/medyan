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
    
    ///Basic monomer consisting of actin, formin, capping, and a virtual front
    class CMonomerBasic : public CMonomer {
        
    private:
        SpeciesFilament* _species[5]; ///<array of contained species
        
    public:
        CMonomerBasic(std::vector<SpeciesFilament*> species, Compartment* c) : CMonomer(c)
        {
            ///Initialize member array of species
            for(auto &s : species) {
                std::string name = s->getName();
                
                ///Set up array
                if(name == "Actin")
                    _species[0] = s;
                
                else if(name == "Front")
                    _species[1] = s;
                
                else if(name == "Back")
                    _species[2] = s;
                
                else if(name == "X-Formin")
                    _species[3] = s;
                
                else if(name == "Capping")
                    _species[4] = s;
                
                else {}
            }
        }
        
        ///Look up a species given a name
        virtual SpeciesFilament* getActin() {return _species[0];}
        virtual SpeciesFilament* getFront() {return _species[1];}
        virtual SpeciesFilament* getBack() {return _species[2];}
        virtual SpeciesFilament* getFormin() {return _species[3];}
        virtual SpeciesFilament* getCapping() {return _species[4];}
        
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
        
    };
    
}





#endif /* defined(__CytoSim__CFilamentImpl__) */
