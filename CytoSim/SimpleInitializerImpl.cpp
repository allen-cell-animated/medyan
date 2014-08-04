//
//  SimpleInitializerImpl.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "SimpleInitializerImpl.h"

namespace chem {
    
    ///Initialize the compartment grid (3D in this implementation)
    void SimpleInitializerImpl::initializeGrid(CompartmentGrid *grid) {
        
        CompartmentGrid temp{NGRID, NGRID, NGRID};
        grid = &temp;
        
        Compartment& cProto = grid->getProtoCompartment();
        
        ///Add species
        cProto.setDiffusionRate(cProto.addSpecies("Actin", 10),_diffusion_rate);
        
        ///Set side length
        std::vector<float> sides{100.0, 100.0, 100.0};
        cProto.setSides(sides.begin());
        
        ///initialize
        grid->initialize();
        grid->activateAll();
        grid->generateDiffusionReactions();
        
        grid->addChemSimReactions();
    }
    
    ///Initializer, inits a cylinder to have actin, virtual back/front, and formin/capping caps.
    CCylinder* SimpleInitializerImpl::createCCylinder(Compartment* c, CCylinder* lastCCylinder, bool extension)
    {
        
        ///maxlength is same length as mcylinder
        int maxlength = L / monomer_size;
        
        ///Set length
        int length;
        if(extension) length = 1; else length = maxlength;
        
        CCylinder* cylinder = new CCylinder(c);
        
        ///remove front from last ccylinder
        if(lastCCylinder != nullptr && !extension) {
            auto front = static_cast<CMonomerBasic*>(lastCCylinder->getCMonomer(maxlength - 1))->getFront();
            if(front != nullptr) front->setN(0);
        }
        else if(extension) {
            auto front = static_cast<CMonomerBasic*>(lastCCylinder->getCMonomer(maxlength - 1))->getFront();
            front->getRSpecies().down();
        }

        ///Add monomers
        for (int index = 0; index < maxlength; index++)
        {
            ///add to cylinder
            CMonomerBasic* m = new CMonomerBasic(c); cylinder->addCMonomer(m);
            
            //Set initial species
            if(index < length) {
                m->getActin()->setN(1);
            }
            
            if(index == length - 1) {
                m->getFront()->setN(1);
            }
            
            if(index == 0 && lastCCylinder == nullptr) m->getBack()->setN(1);

        }
        Filament* parentFilament = nullptr;
        
        ///Callbacks needed
        auto polyCallback = FilamentPolyCallback(parentFilament);
        auto depolyCallback = FilamentDepolyCallback(parentFilament);
        
        auto extensionCallback = FilamentExtensionCallback(parentFilament);
        auto retractionCallback = FilamentRetractionCallback(parentFilament);
        
        
        //Look up diffusing species in this compartment
        Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
        
        ReactionBase *rPoly, *rDepoly;
        
        ///Loop through all spots in cylinder, add poly reactions
        for (int index = 0; index < maxlength; index++) {
            
            ///Monomer and bounds at current index     
            CMonomerBasic *m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index));
            CMonomerBasic *m2 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index+1));
            
            ///extension callback
            if (index == maxlength - 1){
                rPoly = c->addInternal<Reaction,2,0>({m1->getFront(), actinDiffusing},_k_on_plus);
                boost::signals2::shared_connection_block
                rcb1(rPoly->connect(extensionCallback,false));
            }
            ///Typical case
            else {
                ///Add basic polymerization reactions
                rPoly = c->addInternal<Reaction,2,2>({m1->getFront(), actinDiffusing,
                                                      m2->getActin(), m2->getFront()}, _k_on_plus);
                boost::signals2::shared_connection_block
                rcb1(rPoly->connect(polyCallback,false));
            }
            
            rPoly->setAsPolymerizationReaction();
            cylinder->addReaction(rPoly);
        }
        
        ///loop through all spots in subfilament, add depoly reactions
        for (int index = maxlength - 1; index >= 0; index--) {
            
            ///Monomer and bounds at current index
            
            CMonomerBasic *m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index-1));
            CMonomerBasic *m2 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index));
            
            ///Retraction callback
            if(index == 0) {
                rDepoly = c->addInternal<Reaction,2,0>({m2->getFront(), m2->getActin()},
                                                       _k_off_plus);
                boost::signals2::shared_connection_block
                rcb(rDepoly->connect(retractionCallback,false));
            }
            
            ///Typical case
            else {
                /// add basic depolymerization reactions
                rDepoly = c->addInternal<Reaction,2,2>({m2->getFront(), m2->getActin(),
                                                        m1->getFront(), actinDiffusing}, _k_off_plus);
                boost::signals2::shared_connection_block
                rcb(rDepoly->connect(depolyCallback,false));
            }
            cylinder->addReaction(rDepoly);
        }
        
        ///clean up and return
        return cylinder;
    }
    
    ///Remove a cylinder. in this impl, set the front of the new front cylinder
    void SimpleInitializerImpl::removeCCylinder(CCylinder* cylinder)
    {
        
//        ///Set the new front subfilament to have an end species
//        CSubFilament* frontSubFilament = parentFilament->getFrontCSubFilament();
//        
//        ///get end species (if there is one)
//        std::string endName;
//        auto endSpecies = frontSubFilament->backCMonomer()->getActiveEndSpecies();
//        
//        ///Assign end name from previous subfilament
//        if(endSpecies != nullptr)
//            endName= endSpecies->getName();
//        else endName = "Front";
//        
//        ///remove front sub filament
//        parentFilament->removeCSubFilament(frontSubFilament);
//        
//        ///Set new end of new front subfilament
//        parentFilament->getFrontCSubFilament()->frontCMonomer()
//        ->getSpeciesByName(endName)->getRSpecies().up();
        
    } 
        
    ///Constructor, initializes species container
    CMonomerBasic::CMonomerBasic(Compartment* c)
    {
        ///Initialize member array of species
        _species.push_back(c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("Actin"), 0, 1));
        _species.push_back(c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("Front"), 0, 1));
        _species.push_back(c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("Back"), 0, 1));
        
    }
    
    ///Look up species by name
    Species* CMonomerBasic::getSpeciesByName(std::string& name)
    {
        for (auto &s : _species) {
            if(name.find(s->getName()) != std::string::npos) return s;
        }
        
        return nullptr;
    }
    
    ///Find active filament species
    ///@note return null if none
    SpeciesFilament* CMonomerBasic::getActiveFilamentSpecies() {
        
        auto itEnd = _species.begin() + 1;
        
        for (auto it = _species.begin(); it != itEnd; it++)
            if((*it)->getN() == 1) return (*it);
        return nullptr;
    }
    
    ///Find active end species
    ///@note return null if none
    SpeciesFilament* CMonomerBasic::getActiveEndSpecies() {
    
        for (auto it = _species.begin() + 1; it != _species.end(); it++)
            if((*it)->getN() == 1) return (*it);
        return nullptr;
    }
    
    ///Print a species in this filament element
    void CMonomerBasic::print()
    {
        for (auto &s : _species)
            if(s->getN() == 1) std::cout << s->getName().at(0);
    }
    
    
    CBoundBasic::CBoundBasic(Compartment* c)
    {
        ///Initialize member array of species
        _species.push_back(c->addSpeciesBound("Empty", 0, 1));
    }
    
    ///Look up species by name
    Species* CBoundBasic::getSpeciesByName(std::string& name)
    {
        for (auto &s : _species)
            if(s->getName() == name) return s;
        return nullptr;
    }
    
    ///Print a species in this filament element
    void CBoundBasic::print()
    {
        for (auto &s : _species)
            if(s->getN() == 1) std::cout << s->getName().at(0);
    }
        

    
}; //end namespace chem