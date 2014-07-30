//
//  SimpleInitializerImpl.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "SimpleInitializerImpl.h"

namespace chem {
    
    ///Initializer, based on the given simulation
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    CCylinder* SimpleInitializerImpl::createCCylinder(Compartment* c, std::vector<std::string> species, int length)
    {
        int maxlength = L / monomer_size;
        
        CCylinder* cylinder = new CCylinder(c);
        
        ///Add monomers
        for (int index = 0; index < maxlength; index++)
        {
            ///Add species we want
            SpeciesFilament *actin, *capping, *formin, *back, *front;
            SpeciesBound *empty, *myosin, *ma;
            
            //polymer species
            actin = c->addSpeciesFilament("Actin", 0, 1);
            capping = c->addSpeciesFilament("Capping", 0, 1);
            formin = c->addSpeciesFilament("X-Formin", 0, 1);
            
            ///front and back
            back = c->addSpeciesFilament("Back", 0, 1);
            front = c->addSpeciesFilament("Front", 0, 1);
            
            ///bound species
            myosin = c->addSpeciesBound("Myosin", 0, 1);
            ma = c->addSpeciesBound("A-MyosinActin", 0, 1);
            
            //empty
            empty = c->addSpeciesBound("Empty", 0, 1);
            
            //Set initial species
            if(index < length) {
                actin->setN(1);
                empty->setN(1);
            }
            if((index == length - 1) && (length != maxlength)) {
                if (std::find(species.begin(), species.end(), "X-Formin") != species.end())
                    formin->setN(1);
                else front->setN(1);
            }
            
            if(index == 0 && (parentFilament->numCSubFilaments() == 1)) back->setN(1);
            
            ///add to subfilament
            cylinder->addCMonomer(new CMonomerBasic({actin, front, back, formin, capping}, c));
            cylinder->addCBound(new CBoundBasic({empty, myosin, ma}, c));
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
        
        ///Loop through all spots in subfilament, add poly reactions
        for (int index = 0; index < maxlength; index++) {
            
            ///Monomer and bounds at current index
            //CBoundBasic *b1 = static_cast<CBoundBasic*>(subf->getCBound(index));
            CBoundBasic *b2 = static_cast<CBoundBasic*>(subf->getCBound(index+1));
            
            CMonomerBasic *m1 = static_cast<CMonomerBasic*>(subf->getCMonomer(index));
            CMonomerBasic *m2 = static_cast<CMonomerBasic*>(subf->getCMonomer(index+1));
            
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
            subf->addReaction(rPoly);
            
            ///add capping polymerization and depolymerization reactions (+)
            rPoly = c->addInternal<Reaction,2,1>({cappingDiffusing, m1->getFront(),
                m1->getCapping()}, _k_capping_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getCapping(), m1->getFront(),
                cappingDiffusing}, _k_capping_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            subf->addReaction(rPoly);
            subf->addReaction(rDepoly);
            
            ///add formin polymerization and depolymerization reactions (+)
            
            rPoly = c->addInternal<Reaction,2,1>({forminDiffusing, m1->getFront(),
                m1->getFormin()}, _k_formin_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getFormin(), m1->getFront(),
                forminDiffusing}, _k_formin_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            subf->addReaction(rPoly);
            subf->addReaction(rDepoly);
            
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
            subf->addReaction(rPoly);
        
        ///loop through all spots in subfilament, add depoly reactions
        for (int index = maxlength - 1; index >= 0; index--) {
            
            ///Monomer and bounds at current index
            //CBoundBasic *b1 = static_cast<CBoundBasic*>(subf->getCBound(index-1));
            CBoundBasic *b2 = static_cast<CBoundBasic*>(subf->getCBound(index));
            
            CMonomerBasic *m1 = static_cast<CMonomerBasic*>(subf->getCMonomer(index-1));
            CMonomerBasic *m2 = static_cast<CMonomerBasic*>(subf->getCMonomer(index));
            
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
            subf->addReaction(rDepoly);
        }
        
        ///clean up and return
        return subf;
    }
    
    
    
}; //end namespace chem