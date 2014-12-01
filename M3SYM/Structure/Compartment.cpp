
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Compartment.h"

#include "Visitor.h"

Compartment& Compartment::operator=(const Compartment &other) {
    
    _species.clear();
    _internal_reactions.clear();
    _diffusion_reactions.clear();
    other.cloneSpecies(this);
    other.cloneReactions(this);
    _diffusion_rates = other._diffusion_rates;
    return *this;
    // Note that _neighbours is not copied, as well as activated
    
    // Should copy cylinders, beads, etc in future... not clear yet
    
}
    
bool Compartment::apply_impl(SpeciesVisitor &v) {
    for(auto &s : _species.species()) {
        v.visit(s.get());
    }
    return true;
}

bool Compartment::apply_impl(ReactionVisitor &v) {
    for(auto &r : _internal_reactions.reactions()) {
        v.visit(r.get());
    }
    return true;
}

vector<ReactionBase*> Compartment::generateDiffusionReactions(Compartment* C)
{
    vector<ReactionBase*> rxns;
    
    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        int diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0) // Based on a convention that diffusing reactions require positive rates
            continue;
    
        if(C->isActivated()) {
            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);
            ReactionBase *R = new Reaction<1,1>({sp_this.get(),sp_neighbour},diff_rate);
            R->setReactionType(ReactionType::DIFFUSION);
            this->addDiffusionReactionUnique(unique_ptr<ReactionBase>(R));
            rxns.push_back(R);
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}


void Compartment::generateAllDiffusionReactions() {
    if(_activated) {
        for (auto &C: _neighbours)
            generateDiffusionReactions(C);
    }
}

bool operator==(const Compartment& a, const Compartment& b) {
    if(a.numberOfSpecies()!=b.numberOfSpecies() or a.numberOfInternalReactions()!=b.numberOfInternalReactions())
        return false;
    
    if(typeid(a)!=typeid(b))
        return false;
    
    bool spec_bool = false;
    auto sit_pair = mismatch(a._species.species().begin(),a._species.species().end(),b._species.species().begin(),
                                  [](const unique_ptr<Species> &A, const unique_ptr<Species> &B)
                                  {return (*A)==(*B); });
    if(sit_pair.first==a._species.species().end())
        spec_bool=true;
    
    
    bool reac_bool = false;
    auto rit_pair = mismatch(a._internal_reactions.reactions().begin(),a._internal_reactions.reactions().end(),b._internal_reactions.reactions().begin(),
                                  [](const unique_ptr<ReactionBase> &A, const unique_ptr<ReactionBase> &B)
                                  {return (*A)==(*B);});
    if(rit_pair.first==a._internal_reactions.reactions().end())
        reac_bool=true;
    
    return spec_bool && reac_bool;
}