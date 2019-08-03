
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CCylinder.h"

#include "Cylinder.h"

#include "CBound.h"
#include "ChemManager.h"

ChemSim* CCylinder::_chemSim = 0;

/// Default constructor, sets compartment and cylinder
CCylinder::CCylinder(Compartment* C, Cylinder* c)
    : _compartment(C), _pCylinder(c) {
    //set size based on parent cylinder
    _size = SysParams::Geometry().cylinderSize[c->getType()] /
            SysParams::Geometry().monomerSize[c->getType()];
}


CCylinder::CCylinder(const CCylinder& rhs, Compartment* c)
    : _compartment(c), _pCylinder(rhs._pCylinder), _size(rhs._size) {
        
    CCylinder* rhsPtr = const_cast<CCylinder*>(&rhs);
        
    //copy all monomers, bounds
    for(auto &m : rhs._monomers)
        _monomers.emplace_back(m->clone(c));
    
    //copy all internal reactions
    for(auto &r: rhs._internalReactions) {
        ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
        
        if(r->getCBound() != nullptr)
            r->getCBound()->setOffReaction(rxnClone);
        
        addInternalReaction(rxnClone);
    }
    //copy all cross-cylinder reactions
    for(auto it = rhs._crossCylinderReactions.begin();
             it != rhs._crossCylinderReactions.end(); it++) {
        
        for(auto &r : it->second) {

            //copy cbound if any
            ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
            
            if(r->getCBound() != nullptr)
                r->getCBound()->setOffReaction(rxnClone);
            
            addCrossCylinderReaction(it->first, rxnClone);
        }
    }
    //Copy reacting cylinders, Clone reactions where this cylinder is involved
    for(auto &ccyl : rhs._reactingCylinders) {

        //clone reactions
        for(auto &r: ccyl->getCrossCylinderReactions()[rhsPtr]) {
            
            //copy cbound if any
            ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
            
            if(r->getCBound() != nullptr) {
                #ifdef CHECKRXN
                if(dynamic_cast<CMotorGhost*>(r->getCBound())) {
                    CMotorGhost* cm = dynamic_cast<CMotorGhost*>(r->getCBound());
                    cout<<"MotorGhost mID "<<cm->getMotorGhost()->getId()
                    <<" with Rxn "<<rxnClone<<endl;
                }
				#endif
                r->getCBound()->setOffReaction(rxnClone);
            }
            
            ccyl->addCrossCylinderReaction(this, rxnClone);
        }
    }
}

void CCylinder::addInternalReaction(ReactionBase* r) {
    
    //add to compartment and chemsim
    _compartment->addInternalReaction(r);
    _chemSim->addReaction(r);
    
    //add to local reaction list
    _internalReactions.insert(r);
    
    //activate reaction
    r->activateReaction();
}


void CCylinder::removeInternalReaction(ReactionBase* r) {
    
    //remove from internal reaction list
    if (_internalReactions.find(r) != _internalReactions.end()) {
//        std::cout<<"passivate removeInternalReaction"<<endl;
        //passivate reaction, removing from dependents
        r->passivateReaction();
        
        //remove from compartment and chemsim
        _chemSim->removeReaction(r);
        _compartment->removeInternalReaction(r);
        
        _internalReactions.erase(r);
    }
}

void CCylinder::addCrossCylinderReaction(CCylinder* other,
                                         ReactionBase* r) {
    
    //add to compartment and chemsim
    _compartment->addInternalReaction(r);
    _chemSim->addReaction(r);

    #ifdef CHECKRXN
    if(r->getCBound() != nullptr) {
        if(dynamic_cast<CMotorGhost*>(r->getCBound())) {
            CMotorGhost* cm = dynamic_cast<CMotorGhost*>(r->getCBound());
            cout<<"MotorGhost mID "<<cm->getMotorGhost()->getId()
                <<" with Rxn "<<r<<" with RNodeNRM "<<r->getRNode()<<endl;
        }
    }
    #endif

    //add to this reaction map
    _crossCylinderReactions[other].insert(r);
    other->addReactingCylinder(this);
    
    //activate reaction
    r->activateReaction();
}

void CCylinder::addReactingCylinder(CCylinder* other) {
    _reactingCylinders.insert(other);
}

void CCylinder:: removeAllInternalReactions() {
    
    auto tempReactions = _internalReactions;
    for (auto &r : tempReactions) removeInternalReaction(r);
}

void CCylinder::removeCrossCylinderReaction(CCylinder* other,
                                            ReactionBase* r) {
    #ifdef CHECKRXN
    cout<<"Removing cross Cylinder Reaction"<<endl;
	#endif
    auto it = _crossCylinderReactions[other].find(r);
    if(it != _crossCylinderReactions[other].end()) {
	    #ifdef CHECKRXN
    	cout<<"found"<<endl;
		#endif
       
        //erase the reaction
        _crossCylinderReactions[other].erase(it);
//        std::cout<<"passivate removeCrossCylinderReaction"<<endl;

        //passivate reaction, removing from dependents
        r->passivateReaction();
        
        //remove from compartment and chemsim
        _chemSim->removeReaction(r);
        _compartment->removeInternalReaction(r);
        
        //if number of reactions in cross-cylinder
        //has dropped to zero, delete it
        if(_crossCylinderReactions[other].empty()) {
            
            _crossCylinderReactions.erase(other);
            
            //also remove from reacting of other ccylinder
            auto it2 =other->_reactingCylinders.find(this);
            
            if(it2 != other->_reactingCylinders.end())
                other->_reactingCylinders.erase(it2);
        }
    }
}

void CCylinder::removeCrossCylinderReactions(CCylinder* other) {
    
    auto tempReactions = _crossCylinderReactions[other];
    
    for(auto &r : tempReactions)
        removeCrossCylinderReaction(other, r);
}

void CCylinder::removeAllCrossCylinderReactions() {
    
    auto tempMap = _crossCylinderReactions;
    
    for(auto it = tempMap.begin(); it != tempMap.end(); it++)
        removeCrossCylinderReactions(it->first);
}

void CCylinder::removeReactingCylinder(CCylinder* other) {
    
    other->removeCrossCylinderReactions(this);
}

void CCylinder::removeAllReactingCylinders() {
    
    auto tempReactingCylinders = _reactingCylinders;
    
    for(auto &cc : tempReactingCylinders)
        cc->removeCrossCylinderReactions(this);
}

CCylinder::~CCylinder() {
    
    //Remove all reactions owned by this ccylinder
    removeAllInternalReactions();
    removeAllCrossCylinderReactions(); 
    
    //remove all reactions involving this ccylinder
    removeAllReactingCylinders();
    
    //Remove all species
    for(auto &m: _monomers) {
        
        for(int i = 0; i < CMonomer::_numFSpecies[_pCylinder->getType()]; i++) {
            SpeciesFilament* s = m->speciesFilament(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        for(int i = 0; i < CMonomer::_numBSpecies[_pCylinder->getType()]; i++) {
            SpeciesBound* s = m->speciesBound(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
    }
}

void CCylinder::passivatefilcrossreactions(){
    
    for (auto it2=_crossCylinderReactions.begin(); it2!=_crossCylinderReactions.end(); ++it2){
        auto mySet = it2->second;
        for (auto it: mySet) {
            if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::SEVERING
               ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
               ||it->getReactionType() ==ReactionType::AGING)
            {it->passivateReaction();}
        
        }}
//    auto tempReactions = _crossCylinderReactions[this];
//    if(this->getCylinder()->isPlusEnd()){
//        for(auto &it : tempReactions){
//            std::cout<<it->getReactionType()<<endl;
//        }
//    }
//    for(auto &it : tempReactions){
//        if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
//           ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
//           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
//           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
//           ||it->getReactionType() ==ReactionType::SEVERING
//           ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
//           ||it->getReactionType() ==ReactionType::AGING)
//        {it->passivateReaction();}
//    }
}

void CCylinder::activatefilcrossreactions(){    
    for (auto it2=_crossCylinderReactions.begin(); it2!=_crossCylinderReactions.end(); ++it2){
        auto mySet = it2->second;
        for (auto it: mySet) {
            if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::SEVERING
               ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
               ||it->getReactionType() ==ReactionType::AGING)
            {it->activateReaction();}
            
        }}}
void CCylinder::passivatefilreactions(){
    for(auto &it: _internalReactions){
        if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::SEVERING
           ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
           ||it->getReactionType() ==ReactionType::AGING)
        {it->passivateReaction();}}}
void CCylinder::activatefilreactions(){
    for(auto &it: _internalReactions){
        if(it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::SEVERING
           ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
           ||it->getReactionType() ==ReactionType::AGING)
        {it->activateReaction();}}}

vector<ReactionBase*> CCylinder::getAllReactions() {
    
    vector<ReactionBase*> reactions;
    
    //get internal rxns
    for(auto r : _internalReactions) reactions.push_back(r);
    
    //get cross cylinder rxns
    for(auto it = _crossCylinderReactions.begin();
            it != _crossCylinderReactions.end(); it++)
        
        for(auto it2 = it->second.begin();
                it2 != it->second.end(); it2++)
            
            reactions.push_back(*it2);
    
    return reactions;
}

void CCylinder::printCCylinder()
{
    cout << "Compartment:" << _compartment << endl;
    
    cout << "Composition of CCylinder: " << endl;
    for (auto &m : _monomers){
        m->print();
        cout << ":";
    }
    cout << endl;
}

bool CCylinder::isConsistent() {

    int index = 0;
    for(auto &m : _monomers) {
        
        if(!m->isConsistent())
            cout << "CMonomer inconsistency is at index "
                 << index << "." << endl;
        
        index++;
    }
    return true;
}

short CCylinder::getType() {
    
    return _pCylinder->getType();
}

