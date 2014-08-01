//
//  SimpleInitializerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__SimpleInitializerImpl__
#define __Cyto__SimpleInitializerImpl__

#include <iostream>

#include "ChemInitializerImpl.h"
#include "CCylinder.h"
#include "CFilamentElement.h"
#include "CompartmentContainer.h"
#include "Mcommon.h"
#include "common.h"


namespace chem {

    ///SimpleInitializer is a concrete implementation of the ChemInitailizerImpl class, which sets up CCylinders
    ///to have simple actin network interactions
    class SimpleInitializerImpl : public ChemInitializerImpl {

    private:
        ///REACTION RATES
        //basic
        float _k_on_plus = 21.0;
        float _k_off_plus = 1.4;
        float _k_off_minus = 1.4;
        
        //capping
        float _k_capping_on_plus = 50.0;
        float _k_capping_off_plus = 0.06;
        
        //formin
        float _k_formin_on_plus = 10.0;
        float _k_formin_off_plus = 1.4;
        float _k_accel_on_plus = 100.0;
        
        //Diffusion rate
        float _diffusion_rate = 2000.0;
        
        
    public:
        ///Initialize the compartment grid
        virtual void initializeGrid(CompartmentGrid *grid);

        ///Initializer
        ///@param length - starting length of the CCylinder initialized
        ///@param species - list of species to initialize in CCylinder
        virtual CCylinder* createCCylinder(Compartment* c, std::vector<std::string> species, int length);
        
        ///Remove a CCylinder
        virtual void removeCCylinder(CCylinder *cylinder);
        
    };

    ///Basic monomer consisting of actin, formin, capping, and a virtual front/back
    class CMonomerBasic : public CMonomer {
        
    public:
        ///Constructor, initializes species container
        CMonomerBasic(Compartment* c);
        
        ///Destructor
        ~CMonomerBasic();
        
        ///Get a species
        virtual SpeciesFilament* getActin() {return _species[0];}
        
        virtual SpeciesFilament* getFront() {return _species[0];}
        
        virtual SpeciesFilament* getBack() {return _species[1];}
        
        virtual SpeciesFilament* getFormin() {return _species[2];}
        
        virtual SpeciesFilament* getCapping() {return _species[3];}
        
        ///Look up species by name
        virtual Species* getSpeciesByName(std::string& name);
        
        
        ///Check if this monomer is valid
        bool checkSpecies(int sum)
        {
            return true;
            //            int currentSum = 0;
            //            for(auto &s : _species)
            //                currentSum += s->getN();
            //            return currentSum = sum;
        }
        
        ///Print a species in this filament element
        virtual void print();
        
        ///Find active filament species
        ///@note return null if none
        virtual SpeciesFilament* getActiveFilamentSpecies();
        
        ///Find active end species
        ///@note return null if none
        virtual SpeciesFilament* getActiveEndSpecies();
    };
    
    ///Basic bound consisting of myosin, myosin-actin, and a virtual empty species
    class CBoundBasic : public CBound {
        
    public:
        CBoundBasic(Compartment* c);
        
        ///Destructor
        ~CBoundBasic();
        
        ///Look up a species given a name
        virtual SpeciesBound* getEmpty() {return _species[0];}
        
        ///Look up species by name
        virtual Species* getSpeciesByName(std::string& name);
        
        ///Check if this monomer is valid
        bool checkSpecies(int sum)
        {
            return true;
            //        int currentSum = 0;
            //        for(auto &s : _species)
            //            currentSum += s->getN();
            //        return currentSum = sum;
        }
        
        ///Print a species in this filament element
        virtual void print();
    };

    
    
}; // end namespace chem


#endif /* defined(__Cyto__SimpleInitializerImpl__) */
