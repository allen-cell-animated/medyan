//
//  SimpleInitializerImpl.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "SimpleInitializerImpl.h"
    
///Initialize the compartment grid
void SimpleInitializerImpl::initializeGrid() {
    
    Compartment& cProto = CompartmentGrid::Instance(compartmentGridKey())->getProtoCompartment();
    
    ///Add bulk species
    //cProto.setDiffusionRate(cProto.addSpecies("Actin", 10),_diffusion_rate);
    
    
    CompartmentGrid::Instance(compartmentGridKey())->addSpeciesBulk("Actin", 1000U);
    
    ///initialize all compartments with species
    for(auto &c : CompartmentGrid::Instance(compartmentGridKey())->children())
    {
        Compartment *C = static_cast<Compartment*>(c.get());
        *C = cProto;
    }
    
    ///activate all compartments for diffusion
    CompartmentGrid::Instance(compartmentGridKey())->activateAll();
    
//    ///Generate all diffusion reactions
//    for(auto &c : CompartmentGrid::Instance(compartmentGridKey())->children())
//    {
//        Compartment *C = static_cast<Compartment*>(c.get());
//        C->generateAllDiffusionReactions();
//    }
    
    CompartmentGrid::Instance(compartmentGridKey())->addChemSimReactions();
}

///Initializer, inits a cylinder to have actin, and virtual back/front.
CCylinder* SimpleInitializerImpl::createCCylinder(Filament *pf, Compartment* c,
                                                  bool extensionFront, bool extensionBack)
{
    CCylinder* cylinder = new CCylinder(c);
    
    ///maxlength is same length as mcylinder
    int maxlength = SystemParameters::Geometry().cylinderSize / SystemParameters::Geometry().monomerSize;
    
    ///add to cylinder
    for (int index = 0; index < maxlength; index++) {
        CMonomerBasic* m = new CMonomerBasic(c);
        cylinder->addCMonomer(m);
    }
    
    ///get last ccylinder
    CCylinder* lastCCylinder;
 
    ///extension of front
    if(extensionFront) {

        auto m2 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(0));
        SpeciesFilament* front;
        
        ///remove front from last ccylinder, add to current
        front = m2->getFront(); front->getRSpecies().up();
        
        auto actin = m2->getActin();
        actin->getRSpecies().up();
    }

    ///extension of back
    else if(extensionBack) {
        auto m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(maxlength - 1));
        
        SpeciesFilament* back;
        ///remove back from last cylinder, add to current
        back = m1->getBack(); back->getRSpecies().up();
        
        auto actin = m1->getActin();
        actin->getRSpecies().up();
    }

    ///Base case, initialization
    else {
        Cylinder* lastCylinder = pf->getLastCylinder();
        
        ///remove front from last ccylinder, if not null
        if(lastCylinder != nullptr) {
            lastCCylinder = lastCylinder->getCCylinder();
            auto front = static_cast<CMonomerBasic*>(lastCCylinder->getCMonomer(maxlength - 1))->getFront();
            front->getRSpecies().down();
        }
        ///if first ccylinder, set back
        else {
            auto back = static_cast<CMonomerBasic*>(cylinder->getCMonomer(0))->getBack();
            back->getRSpecies().up();
            
        }
        ///Set as full
        for(int index = 0; index < maxlength; index++) {
            
            auto m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index));
            m1->getActin()->getRSpecies().up();

            if (index == maxlength - 1)  {
                m1->getFront()->getRSpecies().up();
            }
        }
    }
    
    ///Callbacks needed
    auto polyCallback = new FilamentPolyCallback(pf);
    auto extensionFrontCallback = new FilamentExtensionFrontCallback(pf);
    auto extensionBackCallback = new FilamentExtensionBackCallback(pf);

    //Look up diffusing species in this compartment
    Species* actinBulk = &CompartmentGrid::Instance(compartmentGridKey())->findSpeciesBulkByName("Actin");
    
    ReactionBase *rPolyPlus, *rPolyMinus;
    
    ///Loop through all spots in cylinder, add poly reactions
    for (int index = 0; index < maxlength; index++) {
        
        ///Monomer and bounds at current index
        CMonomerBasic *m0 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index-1));
        CMonomerBasic *m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index));
        CMonomerBasic *m2 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index+1));
        
        ///Plus end polymerization
        if (index == maxlength - 1) {
            rPolyPlus = c->addInternal<Reaction,2,0>({m1->getFront(), actinBulk},_k_on_plus);
            boost::signals2::shared_connection_block
                rcb1(rPolyPlus->connect(*extensionFrontCallback,false));
        }
        else {
            ///Add basic polymerization reactions
            rPolyPlus = c->addInternal<Reaction,2,2>({m1->getFront(), actinBulk,
                m2->getActin(), m2->getFront()}, _k_on_plus);
            boost::signals2::shared_connection_block
            rcb1(rPolyPlus->connect(*polyCallback,false));
        }

        ///Minus end polymerization
        if(index == 0) {
            rPolyMinus = c->addInternal<Reaction,2,0>({m1->getBack(), actinBulk},_k_on_minus);
            boost::signals2::shared_connection_block
                rcb1(rPolyMinus->connect(*extensionBackCallback,false));
        }
        else {
            ///Add basic polymerization reactions
            rPolyMinus = c->addInternal<Reaction,2,2>({m1->getBack(), actinBulk,
                                                m0->getActin(), m0->getBack()}, _k_on_minus);
            boost::signals2::shared_connection_block
                rcb2(rPolyMinus->connect(*polyCallback,false));
            
        }
        
        rPolyPlus->setAsPolymerizationReaction();
        cylinder->addReaction(rPolyPlus);
        cylinder->addReaction(rPolyMinus);
    }
    
//    ///loop through all spots in subfilament, add depoly reactions
//    for (int index = maxlength - 1; index >= 0; index--) {
//        
//        ///Monomer and bounds at current index
//        
//        CMonomerBasic *m1 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index-1));
//        CMonomerBasic *m2 = static_cast<CMonomerBasic*>(cylinder->getCMonomer(index));
//        
//        ///Retraction callback
//        if(index == 0) {
//            rDepoly = c->addInternal<Reaction,2,0>({m2->getFront(), m2->getActin()},
//                                                   _k_off_plus);
//            boost::signals2::shared_connection_block
//            rcb(rDepoly->connect(retractionCallback,false));
//        }
//        
//        ///Typical case
//        else {
//            /// add basic depolymerization reactions
//            rDepoly = c->addInternal<Reaction,2,2>({m2->getFront(), m2->getActin(),
//                                                    m1->getFront(), actinDiffusing}, _k_off_plus);
//            boost::signals2::shared_connection_block
//            rcb(rDepoly->connect(depolyCallback,false));
//        }
//        cylinder->addReaction(rDepoly);
//    }
    
    cylinder->addChemSimReactions();
    cylinder->updateReactions();

    ///clean up and return
    return cylinder;
}

///Remove a cylinder. in this impl, set the front of the new front cylinder
void SimpleInitializerImpl::removeCCylinder(Filament* pf, bool retractionFront, bool retractionBack)
{   
} 
    
///Constructor, initializes species container
CMonomerBasic::CMonomerBasic(Compartment* c)
{
    ///Initialize member array of species
    std::vector<std::string> species = {"Actin", "Front", "Back"};
    
    for(auto &s : species)
        _species.push_back(c->addSpeciesFilament(
            SpeciesNamesDB::Instance()->generateUniqueName(s), 0, 1));
    
}

///Look up species by name
Species* CMonomerBasic::getSpeciesByName(std::string& name)
{
    for (auto &s : _species) {
        if(name.find(s->getName()) != std::string::npos) return s;
    }
    
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
