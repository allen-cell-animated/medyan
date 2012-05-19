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
#include <unordered_set>
#include <cmath>
#include <initializer_list>


#include "common.h"
#include "Species.h"

namespace chem {

class RNode;

class Reaction {
private:
    std::vector<Species*> _species;
    std::vector<Reaction*> _dependents;
    RNode* _rnode;
    float _rate;
    const unsigned char _m;
    bool _is_signaling; ///< If true, indicates a signal may be sent when a single step of this Reaction occurs
public:
    // Constructor    
    Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate, bool is_signaling = false);    
    Reaction (const Reaction &r) = delete; // no copying (including all derived classes)
    Reaction& operator=(Reaction &r) = delete;  // no assignment (including all derived classes)
    // Destructor
    ~Reaction();
//    // Setters 
    void setRate(float rate) {_rate=rate;}
    void setRnode(RNode *rhs) {_rnode=rhs;} 
//    // Accessors 
    float getRate() const {return _rate;}
    RNode* getRnode() const {return _rnode;} 
    int getProductOfReactants ()  {
        int prod = 1;
        for(auto sit = beginReactants(); sit!=endReactants(); ++sit) 
            prod*=(*sit)->getN();
        return prod;
    }
    
    /// Return true if this Species emits signals on copy number change
    bool isSignaling () const {return _is_signaling;}
    
    /// Set the signaling behavior of this Reaction
    /// @param is the ChemSignal which will call the associated Signal (typically initiated by the 
    /// Gillespie-like simulation algorithm)
    void makeSignaling (ChemSignal &sm);
    
    /// Destroy the signal associated with this Reaction
    /// @param is the ChemSignal which manages signals
    /// @note To start signaling again, makeSignaling(...) needs to be called
    void stopSignaling (ChemSignal &sm);

    //Iterators
    vr_iterator beginAffected() {return _dependents.begin();}
    vr_iterator endAffected() {return _dependents.end();}
    vsp_iterator beginReactants() {return _species.begin();}
    vsp_iterator endReactants() {return _species.begin()+_m;}
    vsp_iterator beginProducts() {return _species.begin()+_m;}
    vsp_iterator endProducts() {return _species.end();}

//    //Fire the Reaction
    void makeStep() {
        for(auto sit = beginReactants(); sit!=endReactants(); ++sit) (*sit)->down();
        for(auto sit = beginProducts(); sit!=endProducts(); ++sit) (*sit)->up();
    }
    float computePropensity () {
        return std::accumulate(beginReactants(), endReactants(), 
                               _rate, 
                               [](float prod, Species *s){ 
                                   return prod*=s->getN();
                               } );
    }
    void passivateReaction();
    void activateReaction();
    void printSelf () ; 
    void printDependents() ;
    std::vector<Reaction*> getAffectedReactions();
    void registerNewDependent(Reaction *r);
    void unregisterDependent(Reaction *r);
};

} // end of namespace 
#endif
