
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

#include "Compartment.h"

#include "CubeSlicing.h"
#include "GController.h"
#include "MathFunctions.h"
using namespace mathfunc;
#include "Visitor.h"

#include "GTriangle.h"
#include "Triangle.h"

void Compartment::getSlicedVolumeArea() {
    // The calculation requires the
    //  - The position calculation of triangles
    //  - The area calculation of triangles
    //  - The unit normal vector of triangles
    // ASSUMPTIONS:
    //  - This compartment is a CUBE
    size_t numTriangle = _triangles.size();
    if(numTriangle) {
        double sumArea = 0.0;
        array<double, 3> sumNormal {};
        array<double, 3> sumPos {};
        for(Triangle* t: _triangles) {
            double area = t->getGTriangle()->getArea();
            vectorIncrease(sumNormal, vectorMultiply(t->getGTriangle()->getUnitNormal(), area));
            vectorIncrease(sumPos, vectorMultiply(t->coordinate, area));
            sumArea += area;
        }
        double oneOverSumArea = 1.0 / sumArea;
        vectorExpand(sumNormal, oneOverSumArea);
        vectorExpand(sumPos, oneOverSumArea);

        PlaneCubeSlicingResult res = planeCubeSlice(
            sumPos, sumNormal,
            {{
                _coords[0] - SysParams::Geometry().compartmentSizeX * 0.5,
                _coords[1] - SysParams::Geometry().compartmentSizeY * 0.5,
                _coords[2] - SysParams::Geometry().compartmentSizeZ * 0.5
            }},
            SysParams::Geometry().compartmentSizeX // Since it is a cube
        );

        _partialVolume = res.volumeIn;
        _partialArea = res.areaIn;
    }
}

Compartment& Compartment::operator=(const Compartment &other) {
    
    _species.clear();
    _internal_reactions.clear();
    _diffusion_reactions.clear();
    other.cloneSpecies(this);
    other.cloneReactions(this);
    _diffusion_rates = other._diffusion_rates;
    
    return *this;
    
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

vector<ReactionBase*> Compartment::generateDiffusionReactions(Compartment* C) {
    // The compartment C and "this" must be neighbors of each other, and
    // "this" must be an active compartment.

    vector<ReactionBase*> rxns;
    
    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        float diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0)  continue;
    
        if(C->isActivated()) {
            // Scale the diffusion rate according to the contacting areas
            size_t idxFwd = _neighborIndex.at(C), idxBwd = C->_neighborIndex.at(this);
            double scaleFactor =
                0.5 * (_partialArea[idxFwd] + C->_partialArea[idxBwd]) /
                GController::getCompartmentArea()[idxFwd / 2];

            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);
            ReactionBase *R = new DiffusionReaction({sp_this.get(),sp_neighbour}, diff_rate * scaleFactor);
            this->addDiffusionReaction(R);
            rxns.push_back(R);
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

vector<ReactionBase*> Compartment::generateAllDiffusionReactions() {
    
    vector<ReactionBase*> rxns;
    
    if(_activated) {
        for (auto &C: _neighbours) {
            auto newRxns = generateDiffusionReactions(C);
            rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

void Compartment::removeDiffusionReactions(ChemSim* chem, Compartment* C)
{
    //look for neighbor's diffusion reactions
    vector<ReactionBase*> to_remove;
    
    for(auto &r : C->_diffusion_reactions.reactions()) {
        
        auto rs = r.get()->rspecies()[1];
        if(rs->getSpecies().getParent() == this) {
            
            r->passivateReaction();
            chem->removeReaction(r.get());
            
            to_remove.push_back(r.get());
        }
    }
    
    //remove them
    for(auto &r : to_remove)
        C->_diffusion_reactions.removeReaction(r);
    
}

void Compartment::removeAllDiffusionReactions(ChemSim* chem) {
    
    //remove all diffusion reactions that this has ownership of
    for(auto &r : _diffusion_reactions.reactions()) {
        r->passivateReaction();
        chem->removeReaction(r.get());
    }
    
    _diffusion_reactions.clear();
    
    //remove neighboring diffusing reactions with this compartment
    for (auto &C: _neighbours)
        removeDiffusionReactions(chem, C);
}


void Compartment::transferSpecies() {
    
    //get active neighbors
    vector<Compartment*> activeNeighbors;
    
    for(auto &neighbor : _neighbours)
        if(neighbor->isActivated())
            activeNeighbors.push_back(neighbor);
    
    assert(activeNeighbors.size() != 0
           && "Cannot transfer species to another compartment... no neighbors are active");
    
    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;
    
    for(auto &sp : _species.species()) {
        
        int copyNumber = sp->getN();
        auto nit = activeNeighbors.begin();
        
        while(copyNumber > 0) {
            sp->down();
            
            //choose a random active neighbor
            auto neighbor = *nit;
            sp_neighbor = neighbor->findSpeciesByName(sp->getName());
            
            //add to list if not already
            auto spit = find(sp_neighbors.begin(),
                             sp_neighbors.end(),
                             sp_neighbor);
            
            if(spit == sp_neighbors.end())
                sp_neighbors.push_back(sp_neighbor);
            
            //increase copy number
            sp_neighbor->up();
            
            //reset if we've looped through
            if(++nit == activeNeighbors.end())
                nit = activeNeighbors.begin();
            copyNumber--;
        }
        
        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
    }
}

void Compartment::activate(ChemSim* chem) {
    
    assert(!_activated && "Compartment is already activated.");
    
    //set marker
    _activated = true;
    
    //add all diffusion reactions
    auto rxns = generateAllDiffusionReactions();
    for(auto &r : rxns) chem->addReaction(r);

}

void Compartment::deactivate(ChemSim* chem) {
    
    //assert no cylinders in this compartment
    assert((_cylinders.size() == 0)
           && "Compartment cannot be deactivated when containing active cylinders.");
    
    assert(_activated && "Compartment is already deactivated.");
    
    //set marker
    _activated = false;
    
    transferSpecies();
    removeAllDiffusionReactions(chem);
}

bool operator==(const Compartment& a, const Compartment& b) {
    if(a.numberOfSpecies()!=b.numberOfSpecies() or
       a.numberOfInternalReactions()!=b.numberOfInternalReactions())
        return false;
    
    if(typeid(a)!=typeid(b))
        return false;
    
    bool spec_bool = false;
    auto sit_pair = mismatch(a._species.species().begin(),
                             a._species.species().end(),
                             b._species.species().begin(),
            [](const unique_ptr<Species> &A, const unique_ptr<Species> &B)
            {return (*A)==(*B); });
    if(sit_pair.first==a._species.species().end())
        spec_bool=true;
    
    
    bool reac_bool = false;
    auto rit_pair = mismatch(a._internal_reactions.reactions().begin(),
                             a._internal_reactions.reactions().end(),
                             b._internal_reactions.reactions().begin(),
            [](const unique_ptr<ReactionBase> &A, const unique_ptr<ReactionBase> &B)
            {return (*A)==(*B);});
    if(rit_pair.first==a._internal_reactions.reactions().end())
        reac_bool=true;
    
    return spec_bool && reac_bool;
}
