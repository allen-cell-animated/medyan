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
#include "ChemCallbacks.h"
#include "common.h"

///SimpleInitializer is a concrete implementation of the ChemInitailizerImpl class, which sets up CCylinders
///to have simple actin network interactions
class SimpleInitializerImpl : public ChemInitializerImpl {

private:
    ///REACTION RATES
    //basic
    float _k_on_plus = 21.0;
    float _k_on_minus = 5.0;
    float _k_off_plus = 1.4;
    float _k_off_minus = 1.4;
    
    //Diffusion rate
    float _diffusion_rate = 2000.0;
    
    
public:
    ///Initialize the compartment grid
    virtual void initializeGrid();

    ///Initializer
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c, bool extensionFront, bool extensionBack);
    
    ///Remove a CCylinder
    virtual void removeCCylinder(Filament* pf, bool retractionFront, bool retractionBack);
    
};

///Basic monomer consisting of actin, formin, capping, and a virtual front/back
class CMonomerBasic : public CMonomer {
    
public:
    ///Constructor, initializes species container
    CMonomerBasic(Compartment* c);
    
    ///Copy constructor, calls base class
    CMonomerBasic(const CMonomerBasic &rhs, Compartment* c) : CMonomer(rhs, c) {};
    
    ///Destructor
    ~CMonomerBasic() {};
    
    /// Move constructor
    CMonomerBasic (CMonomerBasic &&rhs) noexcept : CMonomer(std::move(rhs)) {
    }
    
    /// Regular Assignment
    CMonomerBasic& operator=(const CMonomerBasic& rhs) = delete;
    
    /// Move assignment
    CMonomerBasic& operator=(CMonomerBasic&& rhs)
    {
        CMonomer::operator=(std::move(rhs));
        return *this;
    }
    
    virtual CMonomerBasic* clone(Compartment* c) {
        return new CMonomerBasic(*this, c);
    }
    
    ///Get a species
    virtual SpeciesFilament* getActin() {return _species[0];}
    
    virtual SpeciesFilament* getFront() {return _species[1];}
    
    virtual SpeciesFilament* getBack() {return _species[2];}
    
    
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
    ///Constructor, initializes species container
    CBoundBasic(Compartment* c);
    
    ///Copy constructor, calls base class
    CBoundBasic(const CBoundBasic &rhs, Compartment* c) : CBound(rhs, c) {};
    
    ///Destructor
    ~CBoundBasic() {};
    
    
    /// Move constructor
    CBoundBasic (CBoundBasic &&rhs) noexcept : CBound(std::move(rhs)) {
    }
    
    /// Regular Assignment
    CBoundBasic& operator=(const CBoundBasic& rhs) = delete;
    
    /// Move assignment
    CBoundBasic& operator=(CBoundBasic&& rhs)
    {
        CBound::operator=(std::move(rhs));
        return *this;
    }
    
    virtual CBoundBasic* clone(Compartment* c) {
        return new CBoundBasic(*this, c);
    }

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

#endif /* defined(__Cyto__SimpleInitializerImpl__) */
