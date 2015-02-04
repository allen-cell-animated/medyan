
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Species.h"

#include "Reaction.h"
#include "Composite.h"

SpeciesNamesDB* SpeciesNamesDB::_instance = 0;

SpeciesNamesDB* SpeciesNamesDB::Instance() {
    if(_instance==0)
        _instance = new SpeciesNamesDB;
    return _instance;
}

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
