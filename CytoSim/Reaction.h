//
//  Reaction.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Reaction_h
#define CytoSim_Experimenting_Reaction_h

#include <iostream>
#include <algorithm>
#include <utility>
#include <array>
#include <cmath>

#include <boost/flyweight.hpp>

#include "common.h"
#include "Species.h"

struct OrigReactRates {
    float orig_frate;
    float orig_brate;
    //Constructor
    OrigReactRates (float ofrate, float obrate) :  orig_frate(ofrate), orig_brate(obrate) {}
    //Comparison operator is needed by the boost::flyweight
    bool operator==(const OrigReactRates& orr) const 
    {
        return orr.orig_frate==this->orig_frate && orr.orig_brate==this->orig_brate;
    }
    //Hashing is needed by the boost::flyweight    
    friend std::size_t hash_value(OrigReactRates const& orr){
        std::size_t seed = 0;
        boost::hash_combine(seed,orr.orig_frate);
        boost::hash_combine(seed,orr.orig_brate);
        return seed;
    }
};

class ReactionBase {
protected:
    float _frate;
    float _brate;
    float _af;
    float _ab;
    flyweight<OrigReactRates> _orig_rates;
public:
    // Constructors
    ReactionBase (float frate, float brate) : _frate(frate), _brate(brate), _af(0), _ab(0), _orig_rates(frate,brate)  {}
    ReactionBase (const ReactionBase &r) = delete; // no copying (including all derived classes)
    ReactionBase& operator=(ReactionBase &r) = delete;  // no assignment (including all derived classes)
    // Setters 
    void setFRate(float rate){_frate=rate;}
    void setBRate(float rate){_brate=rate;}
    // Accessors 
    float getFRate(){return _frate;}
    float getBRate(){return _brate;}
    float getOrigFRate(){return _orig_rates.get().orig_frate;}
    float getOrigBRate(){return _orig_rates.get().orig_brate;}
    float getFPropensity() {return _af;}
    float getBPropensity() {return _ab;}
    //Fire Reactions
    void makeFStep() {this->makeFStepImpl(); }
    void makeBStep() {this->makeBStepImpl(); }
    void updateFPropensity() {this->updateFPropensityImpl();}
    void updateBPropensity() {this->updateBPropensityImpl();}
    void updatePropensities() {updateFPropensity(); updateBPropensity();}
    void printSelf() {
        std::cout << "Reaction, ptr=" << this << "\n";
        std::cout << "OrigRates: " << getOrigFRate() << ", " << getOrigBRate() << "\n";
        std::cout << "CurrRates: " << getFRate() << ", " << getBRate() << "\n";

        this->printSelfImpl();
    }
    // (Implementation: Virtual Functions)
protected:
    virtual void makeFStepImpl() = 0;
    virtual void makeBStepImpl() = 0;
    virtual void updateFPropensityImpl() = 0;
    virtual void updateBPropensityImpl() = 0;
    virtual void printSelfImpl() = 0;
    //Destructor
    virtual ~ReactionBase() {};
};


template <unsigned char M,unsigned char N> 
class Reaction : public ReactionBase {
protected:
    std::array<Species*, M> _LSpecies;
    std::array<Species*, N> _RSpecies;
public:
    Reaction<M,N>(float f_rate, float b_rate, std::array<Species*, M> &LSpecies, std::array<Species*, N> &RSpecies) : 
    ReactionBase(f_rate, b_rate), _LSpecies(LSpecies), _RSpecies(RSpecies) {
        for(auto s : LSpecies)
            s->addFReaction(this);
        for(auto s : RSpecies)
            s->addBReaction(this);
    }
    ~Reaction<M,N>(){
        for(auto s : _LSpecies)
            s->removeFReaction(this);
        for(auto s : _RSpecies)
            s->removeBReaction(this);
    }
    void makeFStepImpl() override {
        for(auto s : _LSpecies)
            s->down();
        for(auto s : _RSpecies)
            s->up();

    }
    void makeBStepImpl() override {
        for(auto s : _LSpecies)
            s->up();
        for(auto s : _RSpecies)
            s->down();
    }
    void printSelfImpl() override {
        std::cout << "\nReactants\n";
        for(auto s : _LSpecies)
            s->printSelf();
        std::cout << "Products:\n";
        for(auto s : _RSpecies)
            s->printSelf();
    }
};

template<> class Reaction<1,1> : public ReactionBase {
protected:
    Species& _lspecies;
    Species& _rspecies;
public:
    Reaction<1,1>(float f_rate, float b_rate, Species &lspecies, Species &rspecies) : 
    ReactionBase(f_rate, b_rate), _lspecies(lspecies), _rspecies(rspecies) {
        _lspecies.addFReaction(this);
        _rspecies.addFReaction(this);
    }
    ~Reaction<1,1>(){
        _lspecies.removeFReaction(this);
        _rspecies.removeFReaction(this);
    }
    void updateFPropensityImpl() override {
        _af=_frate*_lspecies.getN();
    }
    void updateBPropensityImpl() override {
        _ab=_brate*_rspecies.getN();
    }
    void makeFStepImpl() override { 
        _lspecies.down();
        _rspecies.up();
        updatePropensities();
    }
    void makeBStepImpl() override { 
        _lspecies.up();
        _rspecies.down();
        updatePropensities();
    }
    void printSelfImpl() override {
        std::cout << "template<> Reaction<1,1>...\n";
        std::cout << "\nReactants\n";
        _lspecies.printSelf();
        std::cout << "Products:\n";
        _rspecies.printSelf();
    }
};

#endif
