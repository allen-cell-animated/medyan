
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Species.h"

#include "Reaction.h"
#include "Composite.h"

#ifdef RSPECIES_SIGNALING
boost::signals2::connection Species::connect(
    function<void (RSpecies *, int)> const &RSpecies_callback, int priority)
{
    if (!isSignaling())
        startSignaling(); 
    return _rspecies->_signal->connect(priority, RSpecies_callback);
}
#endif

Composite* Species::getRoot()
{
    if(hasParent())
       return this->getParent()->getRoot();
    return nullptr;
}

ostream& operator<<(ostream& os, const Species& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

void Species::updateReactantPropensities() {
    
    for(auto r : _rspecies->reactantReactions()){
        r->updatePropensity();
    }
}

void Species::activateReactantReactions() {
    
    for(auto r : _rspecies->reactantReactions())
        r->activateReaction();
}


void Species::passivateReactantReactions() {
    
    for(auto r : _rspecies->reactantReactions())
        r->passivateReaction();
}


unordered_map<string,int> SpeciesNamesDB::_map_string_int;
vector<string> SpeciesNamesDB::_vec_int_string;
unsigned long  SpeciesNamesDB::_num = 0;
