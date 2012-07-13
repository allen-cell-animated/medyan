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


/// Reaction class represents simple chemical reactions of the form A + B -> C.  

/*! Reaction is defined in terms of reactant and product RSpecies. In the current implementation, 
 *  the stoichiometric coefficents can only be one: i.e. A + B -> C + D is allowed, but A + 2B -> C + D is not allowed. 
 *  Almost all chemical reactions are either unimolecular or bimolecular, so this restriction should not be too
 *  burdensom. Also, the Reaction indicates a forward process only. For processes in both directions, e.g. A <-> B, 
 *  two Reactions objects need to be defined, corresponding to A->B and B->A. 
 *
 *  A Reaction tracks other Reaction objects that are affected if this Reaction is executed. A Reaction may be set 
 *  up such that it "signals" when a reaction event happens, in which case the corresponding callbacks are called.
 *
 */
class Reaction {
private:
    std::vector<RSpecies*> _rspecies; ///< Reactants and products constituting this Reaction
    std::vector<Reaction*> _dependents; ///< Pointers to Reaction objects that depend on this Reaction being executed
    RNode* _rnode; ///< A pointer to an RNode object which is used to implement a Gillespie-like algorithm (e.g. NRM)
    float _rate; ///< the rate for this Reaction
    const unsigned char _m; ///< indicates the number of reactants
    bool _is_signaling; ///< If true, indicates a signal may be sent when a single step of this Reaction occurs
    
public:
//    /// Default Constructor produces a Reaction with no reactants or products, zero rate, etc.
//    Reaction () :  _rnode(nullptr), _rate(0.0), _m(0), _is_signaling (false) {}  
//    

    /// The main constructor:
    /// @param Species that are reactants and products are put together into a single list (starting from reactants)
    /// @param M - number of reactants
    /// @param N - number of products
    /// @param rate - the rate constant for this reaction
    Reaction (std::initializer_list<Species*> Species, unsigned char M, unsigned char N, float rate); 
    
    
    Reaction (const Reaction &r) = delete; // no copying (including all derived classes)
    Reaction& operator=(Reaction &r) = delete;  // no assignment (including all derived classes)

    /// Destructor
    ~Reaction();

    /// Sets the reaction rate to the parameter "rate" 
    void setRate(float rate) {_rate=rate;}
    
    /// Sets the RNode pointer associated with this Reaction to rhs. Usually is called only by the 
    /// Gillespie-like algorithms.
    void setRnode(RNode *rhs) {_rnode=rhs;} 

    /// Returns the rate associated with this Reaction
    float getRate() const {return _rate;}
    
    /// Returns a pointer to the RNode associated with this Reaction
    RNode* getRnode() const {return _rnode;} 
    
    /// Returns the number of reactant RSpecies
    unsigned char getM() const {return _m;}
    
    /// Returns the number of product RSpecies
    unsigned char getN() const {return static_cast<unsigned char>(_rspecies.size()-_m);}
    
    /// Computes the product of the copy number of all reactant RSpecies. Can be used to quickly 
    /// determined whether this Reaction is active - i.e. if one of the reactants has a 0 copy number, 
    /// then the propensity for this Reaction is 0.
    int getProductOfReactants () const {
        int prod = 1;
        for(auto sit = cbeginReactants(); sit!=cendReactants(); ++sit) 
            prod*=(*sit)->getN();
        return prod;
    }
    
    /// Return true if this RSpecies emits signals on copy number change
    bool isSignaling () const {return _is_signaling;}
    
    /// Set the signaling behavior of this Reaction
    /// @param is the ChemSignal which will call the associated Signal (typically initiated by the 
    /// Gillespie-like simulation algorithm)
    void makeSignaling (ChemSignal &sm);
    
    /// Destroy the signal associated with this Reaction
    /// @param is the ChemSignal which manages signals
    /// @note To start signaling again, makeSignaling(...) needs to be called
    void stopSignaling (ChemSignal &sm);

    /// Return a const reference to the vector of dependent reactions
    /// @note One can obtain two different lists of affected reactions:
    /// 1) via getAffectedReactions(), where the copy numbers do influence the 
    /// dependencies, and 2) via dependents(), where dependencies stop being counted 
    /// if the copy numbers of reactant species drop to 0. 
    const std::vector<Reaction*>& dependents() {return _dependents;}
        
    //Return an iterator to the beginning of the sequence of reactants
    vrsp_iterator        beginReactants() {return _rspecies.begin();}

    //Return a const iterator to the beginning of the sequence of reactants
    vrsp_const_iterator  cbeginReactants() const {return _rspecies.cbegin();}

    //Return an iterator to the end of the sequence of reactants
    vrsp_iterator        endReactants() {return _rspecies.begin()+_m;}

    //Return a const iterator to the end of the sequence of reactants
    vrsp_const_iterator  cendReactants() const {return _rspecies.cbegin()+_m;}

    //Return an iterator to the beginning of the sequence of products
    vrsp_iterator        beginProducts() {return _rspecies.begin()+_m;}

    //Return a const iterator to the beginning of the sequence of products

    vrsp_const_iterator  cbeginProducts() const {return _rspecies.cbegin()+_m;}

    //Return an iterator to the end of the sequence of products
    vrsp_iterator        endProducts() {return _rspecies.end();}

    //Return a const iterator to the end of the sequence of products
    vrsp_const_iterator  cendProducts() const {return _rspecies.cend();}

    /// Fire the Reaction - make a single step, where reactant RSpecies copy numbers are
    /// decreased by one, and the product RSpecies copy numbers are increased by one.
    /// @note This method does not send a reaction event Signal. The latter is usually done
    /// from within a Gillespie-like algorithm. 
    void makeStep() {
        for(auto sit = beginReactants(); sit!=endReactants(); ++sit) 
            (*sit)->down();
        for(auto sit = beginProducts(); sit!=endProducts(); ++sit) 
            (*sit)->up();
    }
    
    /// Compute the Reaction propensity that is needed by a Gillespie like algorithm:
    /// rate*reactant_1.getN()*reactant_2.getN()...
    float computePropensity () const {
        return std::accumulate(cbeginReactants(), cendReactants(), 
                               _rate, 
                               [](float prod, RSpecies *s){ 
                                   return prod*=s->getN();
                               } );
    }
    
    /// Usually is applied to Reaction objects with propensity of 0 (e.g. when one of the 
    /// copy numbers of reactants has dropped to 0. This method call notifies all other 
    /// Reaction objects that may affect this Reaction to stop tracking this Reaction. Eventually, 
    /// activateReaction() may be called to restart tracking, if the propensity stops being 0.
    void passivateReaction();
    
    /// Requests that Reaction objects that may affect this Reaction to start tracking it, which can be 
    /// used to follow Reaction objects whose propensities change upong firing of some Reaction.
    void activateReaction();
    
    /// Print the Reactions that are affacted by this Reaction being fired
    void printDependents() ;
    
    /// Return the list of Reaction objects that are affected when this Reaction is fired
    /// @return a vector of pointers to the affected Reaction objects
    /// @note This method is "expensive" because it computes from scratch the dependencies. Importantly, 
    /// the copy numbers of molecules do not influence the result of this function. \sa dependents()
    std::vector<Reaction*> getAffectedReactions();
    
    /// Request that the Reaction *r adds this Reaction to its list of dependents which it affects. 
    void registerNewDependent(Reaction *r);
    
    /// Request that the Reaction *r removes this Reaction from its list of dependents which it affects. 
    /// This is usually requested when the reaction propensity drops to zero (i.e. via passivateReaction()).
    void unregisterDependent(Reaction *r);

    /// Print self into an iostream
    friend std::ostream& operator<<(std::ostream& os, const Reaction& rr);    
};

} // end of namespace 
#endif
