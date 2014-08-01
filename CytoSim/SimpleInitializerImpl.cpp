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
        cProto.setDiffusionRate(cProto.addSpecies("Capping", 5),_diffusion_rate);
        cProto.setDiffusionRate(cProto.addSpecies("X-Formin", 5),_diffusion_rate);
        
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
    CCylinder* SimpleInitializerImpl::createCCylinder(Compartment* c, std::vector<std::string> species, int length)
    {
        int maxlength = L / monomer_size;
        
        CCylinder* cylinder = new CCylinder(c);
        
        ///Add monomers
        for (int index = 0; index < maxlength; index++)
        {
            ///add to cylinder
            CMonomerBasic* m = new CMonomerBasic(c); cylinder->addCMonomer(m);
            CBoundBasic* b = new CBoundBasic(c); cylinder->addCBound(b);
            
            //Set initial species
            if(index < length) {
                m->getActin()->setN(1);
                b->getEmpty()->setN(1);
            }
            
            if((index == length - 1) && (length != maxlength)) {
                if (std::find(species.begin(), species.end(), "X-Formin") != species.end())
                    m->getFormin()->setN(1);
                else m->getFront()->setN(1);
            }
            
            if(index == 0 && (parentFilament->numCSubFilaments() == 1)) m->getBack()->setN(1);

        }
        
//        ///Callbacks needed
//        auto polyCallback =
//        CFilamentPolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem, parentFilament);
//        
//        auto depolyCallback =
//        CFilamentDepolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem, parentFilament);
//        
//        auto extensionCallback =
//        CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
//                                         parentFilament, {"Actin"});
//        auto extensionCallback2 =
//        CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
//                                         parentFilament, {"Actin", "Formin"});
//        auto retractionCallback =
//        CFilamentRetractionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem, parentFilament);
        
        
        //Look up diffusing species in this compartment
        Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
        Species* cappingDiffusing = c->findSpeciesDiffusingByName("Capping");
        Species* forminDiffusing = c->findSpeciesDiffusingByName("X-Formin");
        
        ReactionBase *rPoly, *rDepoly;
        
        ///Loop through all spots in cylinder, add poly reactions
        for (int index = 0; index < maxlength; index++) {
            
            ///Monomer and bounds at current index
            CBoundBasic *b2 = static_cast<CBoundBasic*>(cylinder->getCBound(index+1));
            
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
                rPoly = c->addInternal<Reaction,2,3>({m1->getFront(), actinDiffusing, m2->getActin(),
                    m2->getFront(), b2->getEmpty()}, _k_on_plus);
                boost::signals2::shared_connection_block
                rcb1(rPoly->connect(polyCallback,false));
            }
            
            rPoly->setAsPolymerizationReaction();
            cylinder->addReaction(rPoly);
            
            ///add capping polymerization and depolymerization reactions (+)
            rPoly = c->addInternal<Reaction,2,1>({cappingDiffusing, m1->getFront(),
                m1->getCapping()}, _k_capping_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getCapping(), m1->getFront(),
                cappingDiffusing}, _k_capping_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            cylinder->addReaction(rPoly);
            cylinder->addReaction(rDepoly);
            
            ///add formin polymerization and depolymerization reactions (+)
            
            rPoly = c->addInternal<Reaction,2,1>({forminDiffusing, m1->getFront(),
                m1->getFormin()}, _k_formin_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getFormin(), m1->getFront(),
                forminDiffusing}, _k_formin_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            cylinder->addReaction(rPoly);
            cylinder->addReaction(rDepoly);
            
            ///add accelerated polymerization of actin with anti-cappped end
            if (index < maxlength - 1) {
                rPoly =
                c->addInternal<Reaction,2,3>({actinDiffusing, m1->getFormin(), m2->getActin(),
                    m2->getFormin(), b2->getEmpty()}, _k_accel_on_plus);
                
                boost::signals2::shared_connection_block
                rcb(rPoly->connect(polyCallback, false));
            }
            ///extension callback
            else {
                rPoly = c->addInternal<Reaction,2,0>({m1->getFormin(), actinDiffusing},
                                                     _k_accel_on_plus);
                
                boost::signals2::shared_connection_block
                rcb(rPoly->connect(extensionCallback2,false));
            }
            
            rPoly->setAsPolymerizationReaction();
            cylinder->addReaction(rPoly);
        
        }
        ///loop through all spots in subfilament, add depoly reactions
        for (int index = maxlength - 1; index >= 0; index--) {
            
            ///Monomer and bounds at current index
            //CBoundBasic *b1 = static_cast<CBoundBasic*>(subf->getCBound(index-1));
            CBoundBasic *b2 = static_cast<CBoundBasic*>(cylinder->getCBound(index));
            
            CMonomerBasic *m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index-1));
            CMonomerBasic *m2 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index));
            
            ///Retraction callback
            if(index == 0) {
                rDepoly = c->addInternal<Reaction,3,0>({m2->getFront(), m2->getActin(), b2->getEmpty()},
                                                       _k_off_plus);
                boost::signals2::shared_connection_block
                rcb(rDepoly->connect(retractionCallback,false));
            }
            
            ///Typical case
            else {
                /// add basic depolymerization reactions
                rDepoly = c->addInternal<Reaction,3,2>({m2->getFront(), m2->getActin(), b2->getEmpty(),
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
        _species.push_back(c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("Capping"), 0, 1));
        _species.push_back(c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("Formin"), 0, 1));
        
    }
    
    ///Look up species by name
    Species* CMonomerBasic::getSpeciesByName(std::string& name)
    {
        for (auto &s : _species)
            if(s->getName() == name) return s;
        
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