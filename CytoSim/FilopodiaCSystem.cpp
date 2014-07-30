//
//  FilopodiaCSystem.cpp
//  CytoSim
//
//  Created by James Komianos on 7/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CMembrane.h"

namespace chem {
    
    ///Find the current polymerization reactions associated with this CFilament
    template<size_t NDIM>
    std::vector<ReactionBase*>
        SimpleInitializer<NDIM>::findPolymerizationReactions(CFilament* f)
    {
        std::vector<ReactionBase*> polyReactions;
        
        ///go to the front subfilament, front monomer;
        CMonomer* frontCMonomer = f->getLeadingCMonomer();
        
        ///find all reactions associated with this active end polymerization
        Species* end = frontCMonomer->getActiveEndSpecies();
        auto reactions = end->getRSpecies().ReactantReactions();
        
        for(auto it = reactions.begin(); it != reactions.end(); it++)
            if((*it)->isPolymerizationReaction())
                polyReactions.push_back(*it);
        
        return std::vector<ReactionBase*>(polyReactions.begin(), polyReactions.end());
    };
    
    
    ///Update filaments based on a reaction
    ///In this implementation, update polymerization rates based on membrane
    template<size_t NDIM>
    void SimpleInitializer<NDIM>::update(CFilament* f, ReactionBase* r)
    {
        
        //update associated filament polymerization reactions
        _membrane.updateFilamentReactions(f,
                SimpleInitializer<NDIM>::findPolymerizationReactions(f));
        
        //recalculate rates
        _membrane.updateHeight();
        _membrane.updateRates();
    }

    
    ///Initialize proto compartment based on this implementation
    template <size_t NDIM>
    void SimpleInitializer<NDIM>::initializeProtoCompartment(CompartmentSpatial<NDIM>& Cproto)
    {
        ///Add species
        Cproto.setDiffusionRate(Cproto.addSpecies("Actin", 10),_diffusion_rate);
        Cproto.setDiffusionRate(Cproto.addSpecies("Capping", 5),_diffusion_rate);
        Cproto.setDiffusionRate(Cproto.addSpecies("X-Formin", 5),_diffusion_rate);
        Cproto.setDiffusionRate(Cproto.addSpecies("Myosin", 0),_diffusion_rate);
        
        ///Set side length
        std::vector<float> sides{_side_length};
        Cproto.setSides(sides.begin());
    }
    
    ///Remove a CSubFilament, based on the given simulation
    template <size_t NDIM>
    void SimpleInitializer<NDIM>::removeCSubFilament(CFilament* parentFilament)
    {
        ///Set the new front subfilament to have an end species
        CSubFilament* frontSubFilament = parentFilament->getFrontCSubFilament();

        ///get end species (if there is one)
        std::string endName;
        auto endSpecies = frontSubFilament->backCMonomer()->getActiveEndSpecies();
        
        ///Assign end name from previous subfilament
        if(endSpecies != nullptr)
            endName= endSpecies->getName();
        else endName = "Front";
        
        ///remove front sub filament
        parentFilament->removeCSubFilament(frontSubFilament);
        
        ///Set new end of new front subfilament
        parentFilament->getFrontCSubFilament()->frontCMonomer()
                    ->getSpeciesByName(endName)->getRSpecies().up();
    }
    
    //Initialize a number of filaments
    template <size_t NDIM>
    CFilament* FilopodiaCSystem<NDIM>::initializeCFilament(int length)
    {
        CFilament* f = new CFilament();
        CompartmentSpatial<NDIM>* cStart;
        
        ///Starting compartment for 1D, all filaments start in compartment 0
        if (NDIM == 1) {
            cStart = CSystem<NDIM>::_grid->getCompartment(0);
        }
        else {
            std::cout << "Multiple dimensional implementation \
                    not optional for Filopodia (yet). Exiting." << std::endl;
            exit(EXIT_FAILURE);
        }
        
        ///number of units to initialize in this filament
        int numUnits = length / monomer_size;
        int maxUnits = f->maxLengthSubFilament();
        
        ///initialize each filament
        Compartment* cNext = cStart;
        
        int numSubFilaments = (numUnits - 1) / maxUnits + 1;
        int ci = 1;
        
        for(int si = 0; si < numSubFilaments; si++) {
            
            ///Add sub filament to compartment
            cNext->addSubFilament();
            
            int setLength; ///length to intialize subfilament to
            
            if (si == numSubFilaments - 1) 
                setLength = numUnits - (maxUnits * (numSubFilaments - 1));
            else
                setLength = maxUnits;
            
            ///create a subfilament
            CSystem<NDIM>::_initializer->createCSubFilament(f, cNext, {"Actin"}, setLength);
            
            ///move to next compartment if needed
            if(float(si * maxUnits) * monomer_size >= cStart->sides()[0] * float(ci)) {
                
                cNext = cNext->neighbours().back();
                f->increaseNumCompartments();
                ci++;
            }
        }
        f->setLength(numUnits);
        CSystem<NDIM>::_filaments.emplace_back(f);
        CSystem<NDIM>::update(f, nullptr);
            
        return CSystem<NDIM>::_filaments.back().get();
    }

    ///Extend the front of a filament
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::extendFrontOfCFilament(CFilament *f,
                                                        std::vector<std::string> species)
    {
        ///Find next compartment (1D for now)
        Compartment* cCurrent = f->getFrontCSubFilament()->compartment();
        Compartment* cNext;
        
        float compartmentLength = CSystem<NDIM>::_grid->getProtoCompartment().sides()[0];
        
        if(f->length() * monomer_size >= f->numCompartments() * compartmentLength) {
            
            if(cCurrent->neighbours().size() == 1 && cCurrent != CSystem<NDIM>::_grid->getCompartment(0)) {
                std::cout << "A filament has reached the end of the grid. Exiting." << std::endl;
                std::cout << std::endl;
                CSystem<NDIM>::printFilaments();
                std::cout << "tau=" << tau() <<std::endl;
                std::cout << "Done!" <<std::endl;
                exit(EXIT_SUCCESS);
            }
            
            cNext = cCurrent->neighbours().back();
            f->increaseNumCompartments();
            
            ///activate the new compartment if needed
            if(!cNext->isActivated()) {
                
                cNext->addSubFilament();
                
                ///generate diffusion reactions
                for(auto &neighbor : cNext->neighbours()) {
                    
                    auto diffusionReactions = cNext->generateDiffusionReactions(neighbor);
                    
                    for(auto &r : diffusionReactions) {
                        CSystem<NDIM>::_initializer->getChemSim().addReaction(r);
                        r->activateReaction();
                    }
                    
                    diffusionReactions = neighbor->generateDiffusionReactions(cNext);
                    
                    for(auto &r : diffusionReactions) {
                        CSystem<NDIM>::_initializer->getChemSim().addReaction(r);
                        r->activateReaction();
                    }
                }
                
//                ///move half of diffusing species to the new compartment
//                for(auto &sp_this : cNext->species.species()) {
//                    int molecule = sp_this->getMolecule();
//                    int diff_rate = _diffusion_rates[molecule];
//                    if(diff_rate<0) // Based on a convention that diffusing reactions require positive rates
//                        continue;
//                }
            }   
        }
        
        else {
            cCurrent->addSubFilament();
            cNext = cCurrent;
        }
        
        ///Initialize new subfilament
        CSubFilament* s = CSystem<NDIM>::_initializer->createCSubFilament(f, cNext, species, 1);
        
        ///Increase length
        f->increaseLength();
        
        ///Add all new reactions
        for(auto &r : s->getReactions())
            CSystem<NDIM>::_initializer->getChemSim().addReaction(r);
        
    }
    
    ///Retract the front of a CFilament
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::retractFrontOfCFilament(CFilament *f)
    {
        ///remove subfilament from compartment
        CSubFilament* s = f->getFrontCSubFilament();
        Compartment* c = s->compartment();

        c->removeSubFilament();
    
        ///if num sub filaments is one, remove filament from system
        if(f->numCSubFilaments() == 1) {
           
            return;
            
//            auto child_iter = std::find_if(CSystem<NDIM>::_filaments.begin(),CSystem<NDIM>::_filaments.end(),
//                                           [f](const std::unique_ptr<CFilament> &element)
//                                           {
//                                               return element.get()==f ? true : false;
//                                           });
//            if(child_iter!=CSystem<NDIM>::_filaments.end())
//                CSystem<NDIM>::_filaments.erase(child_iter);
//            return;
        }
        
        ///Remove subfilament
        CSystem<NDIM>::_initializer->removeCSubFilament(f);
        
        ///Decrease length
        f->decreaseLength();
    }
  
    
    ///perform one step of retrograde flow
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::retrogradeFlow()
    {
        
        ///loop through all filaments, push species back
        for(auto it = CSystem<NDIM>::_filaments.begin();
                 it != CSystem<NDIM>::_filaments.end(); it++) {
            
            CFilament* f = (*it).get();
            CSubFilament* s = f->getFrontCSubFilament();
            
            ///retract if needed
            if (f->lengthFrontSubFilament() == 1)
                retractFrontOfCFilament(f);
            
            else {
                int currentIndex = f->lengthFrontSubFilament() - 1;
                
                ///get end species name
                Species* endSpecies = s->getCMonomer(currentIndex)->getActiveEndSpecies();
                std::string endName = endSpecies->getName();
                
                //decrease copy numbers
                endSpecies->getRSpecies().down();
                s->getCMonomer(currentIndex)->
                    getActiveFilamentSpecies()->getRSpecies().down();
                
                ///set copy number of new end
                s->getCMonomer(currentIndex - 1)->
                    getSpeciesByName(endName)->getRSpecies().up();
                
                ///decrease length
                f->decreaseLength();
            }
            
        }
    }
    
    //Specializations
    template std::vector<ReactionBase*> SimpleInitializer<1>::findPolymerizationReactions(CFilament* f);
    template std::vector<ReactionBase*> SimpleInitializer<2>::findPolymerizationReactions(CFilament* f);
    template std::vector<ReactionBase*> SimpleInitializer<3>::findPolymerizationReactions(CFilament* f);
    
    template void SimpleInitializer<1>::initializeProtoCompartment(CompartmentSpatial<1>& Cproto);
    template void SimpleInitializer<2>::initializeProtoCompartment(CompartmentSpatial<2>& Cproto);
    template void SimpleInitializer<3>::initializeProtoCompartment(CompartmentSpatial<3>& Cproto);
    
    template void SimpleInitializer<1>::update(CFilament* f, ReactionBase* r);
    template void SimpleInitializer<2>::update(CFilament* f, ReactionBase* r);
    template void SimpleInitializer<3>::update(CFilament* f, ReactionBase* r);
    
    template CSubFilament*
    SimpleInitializer<1>::createCSubFilament(CFilament* parentFilament, Compartment* c, std::vector<std::string> species, int length);
    template CSubFilament*
    SimpleInitializer<2>::createCSubFilament(CFilament* parentFilament, Compartment* c, std::vector<std::string> species, int length);
    template CSubFilament*
    SimpleInitializer<3>::createCSubFilament(CFilament* parentFilament, Compartment* c, std::vector<std::string> species, int length);
    
    
    template void SimpleInitializer<1>::removeCSubFilament(CFilament* parentFilament);
    template void SimpleInitializer<2>::removeCSubFilament(CFilament* parentFilament);
    template void SimpleInitializer<3>::removeCSubFilament(CFilament* parentFilament);
    
    
    //Specializations
    template CFilament* FilopodiaCSystem<1>::initializeCFilament(int length);
    template CFilament* FilopodiaCSystem<2>::initializeCFilament(int length);
    template CFilament* FilopodiaCSystem<3>::initializeCFilament(int length);
    
    template void FilopodiaCSystem<1>::extendFrontOfCFilament(CFilament *f, std::vector<std::string> species);
    template void FilopodiaCSystem<2>::extendFrontOfCFilament(CFilament *f, std::vector<std::string> species);
    template void FilopodiaCSystem<3>::extendFrontOfCFilament(CFilament *f, std::vector<std::string> species);
    
    template void FilopodiaCSystem<1>::retractFrontOfCFilament(CFilament* f);
    template void FilopodiaCSystem<2>::retractFrontOfCFilament(CFilament* f);
    template void FilopodiaCSystem<3>::retractFrontOfCFilament(CFilament* f);
    
    template void FilopodiaCSystem<1>::retrogradeFlow();
    template void FilopodiaCSystem<2>::retrogradeFlow();
    template void FilopodiaCSystem<3>::retrogradeFlow();
    
    
}; //end namespace chem


