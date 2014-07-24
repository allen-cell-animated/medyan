//
//  FilopodiaCSystem.cpp
//  CytoSim
//
//  Created by James Komianos on 7/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "FilopodiaCSystem.h"
#include "CMembrane.h"

namespace chem {
    
    ///Find the current polymerization reactions associated with this CFilament
    template<size_t NDIM>
    std::vector<ReactionBase*>
        SimpleInitializer<NDIM>::findPolymerizationReactions(CFilament* f)
    {
        std::vector<ReactionBase*> polyReactions;
        ///go to the front subfilament, front monomer
        CMonomer* frontCMonomer = f->getFrontCSubFilament()->frontCMonomer();
        
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
        //update associated filament reactions
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
        Cproto.setDiffusionRate(Cproto.addSpecies("Actin",10U),_diffusion_rate);
        Cproto.setDiffusionRate(Cproto.addSpecies("Capping",0U),_diffusion_rate);
        Cproto.setDiffusionRate(Cproto.addSpecies("X-Formin",0U),_diffusion_rate);
        Cproto.setDiffusionRate(Cproto.addSpecies("Myosin", 5U),_diffusion_rate);
        
        ///Set side length
        std::vector<float> sides{_side_length};
        Cproto.setSides(sides.begin());
    }
    
    
    ///Initializer, based on the given simulation
    ///@param length - starting length of the filament initialized
    ///@param maxlength - length of entire reactive filament
    template <size_t NDIM>
    CSubFilament* SimpleInitializer<NDIM>::createCSubFilament(CFilament* parentFilament,
                                                               Compartment* c,
                                                               std::vector<std::string> species,
                                                               int length)
    {
        int maxlength = L / monomer_size;
        
        CSubFilament* subf = new CSubFilament(c);
        parentFilament->addCSubFilament(subf);
        
        ///Add monomers
        for (int index = 0; index < maxlength; index++)
        {
            ///Add species we want
            SpeciesFilament *actin, *capping, *formin, *back, *front;
            SpeciesBound *empty, *myosin, *ma;
            
            //polymer species
            actin = c->addSpeciesFilament("Actin",0,1);
            capping = c->addSpeciesFilament("Capping",0,1);
            formin = c->addSpeciesFilament("X-Formin",0,1);
            
            ///front and back
            back = c->addSpeciesFilament("Back",0,1); front = c->addSpeciesFilament("Front",0,1);
            ///bound species
            myosin = c->addSpeciesBound("Myosin",0,1); ma = c->addSpeciesBound("A-MyosinActin",0,1);
            //empty
            empty = c->addSpeciesBound("Empty",0,1);
            
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
            subf->addCMonomer(new CMonomerBasic({actin, front, back, formin, capping}, c));
            subf->addCBound(new CBoundBasic({empty, myosin, ma}, c));
        }
        
        ///Callbacks needed
        
        auto polyCallback =
            CFilamentPolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                        parentFilament);
        
        auto depolyCallback =
            CFilamentDepolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                          parentFilament);
        
        auto extensionCallback =
            CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                             parentFilament, {"Actin"});
        
        auto extensionCallback2 =
            CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                             parentFilament, {"Actin", "Formin"});
        
        auto retractionCallback =
            CFilamentRetractionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                              parentFilament);
        
        
        //Look up diffusing species in this compartment
        Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
        Species* cappingDiffusing = c->findSpeciesDiffusingByName("Capping");
        Species* forminDiffusing = c->findSpeciesDiffusingByName("X-Formin");
        //Species* myosinDiffusing = c->findSpeciesDiffusingByName("Myosin");
        
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
            CFilamentInitializer<NDIM>::_chem.addReaction(rPoly);
        
            ///add capping polymerization and depolymerization reactions (+)
            rPoly = c->addInternal<Reaction,2,1>({cappingDiffusing, m1->getFront(),
                                                  m1->getCapping()}, _k_capping_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getCapping(), m1->getFront(),
                                                    cappingDiffusing}, _k_capping_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
        
            ///add formin polymerization and depolymerization reactions (+)
            
            rPoly = c->addInternal<Reaction,2,1>({forminDiffusing, m1->getFront(),
                                                  m1->getFormin()}, _k_formin_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getFormin(), m1->getFront(),
                                                    forminDiffusing}, _k_formin_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
            
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
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
        
//            ///add motor binding and unbinding, loading and unloading
//            
//            ReactionBase *rBinding =
//                c->addInternal<Reaction,2,1>({myosinDiffusing, b1->getEmpty(),
//                                              b1->getMyosin()}, _k_binding);
//            
//            ReactionBase *rUnbinding =
//                c->addInternal<Reaction,1,2>({b1->getMyosin(), b1->getEmpty(),
//                                              myosinDiffusing}, _k_unbinding);
//            
//            ReactionBase *rLoading =
//                c->addInternal<Reaction,2,1>({actinDiffusing, b1->getMyosin(),
//                                              b1->getMyosinActin()}, _k_load);
//            
//            ReactionBase *rUnloading =
//                c->addInternal<Reaction,1,2>({b1->getMyosinActin(), b1->getMyosin(),
//                                              actinDiffusing}, _k_unload);
//            
//            
//            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rBinding);
//            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rUnbinding);
//            
//            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rLoading);
//            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rUnloading);
//
//            //add motor stepping
//            ReactionBase *rMForwardStep;
//            ReactionBase *rMBackwardStep;
//            
//            if (index < maxlength - 1) {
//                rMForwardStep =
//                    c->addInternal<Reaction,2,2>({b1->getMyosin(), b2->getEmpty(),
//                                                  b2->getMyosin(), b1->getEmpty()},
//                                                 _k_forward_step);
//                rMBackwardStep =
//                    c->addInternal<Reaction,2,2>({b2->getMyosin(), b1->getEmpty(),
//                                                  b1->getMyosin(), b2->getEmpty()},
//                                                 _k_forward_step);
//                
//                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMBackwardStep);
//                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMForwardStep);
//            }
//            
//            ReactionBase *rMA3ForwardStep;
//            ReactionBase *rMA3BackwardStep;
//            
//            if (index < maxlength - 1) {
//                rMA3ForwardStep =
//                    c->addInternal<Reaction,2,2>({b1->getMyosinActin(), b2->getEmpty(),
//                                                  b2->getMyosinActin(), b1->getEmpty()},
//                                                 _k_forward_step);
//                
//                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3ForwardStep);
//
//                rMA3BackwardStep =
//                    c->addInternal<Reaction,2,2>({b2->getMyosinActin(), b1->getEmpty(),
//                                                  b1->getMyosinActin(), b2->getEmpty()},
//                                                  _k_forward_step);
//                
//                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3BackwardStep);
//            }
//            
//            ReactionBase *rMA3Unbinding =
//                c->addInternal<Reaction,1,3>({b1->getMyosinActin(), myosinDiffusing,
//                                              b1->getEmpty(), actinDiffusing},
//                                              _k_unbinding);
//            
//            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3Unbinding);
        }
        
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
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
            
        }
        ///clean up and return
        return subf;
    }

    //Initialize a number of filaments
    template <size_t NDIM>
    CFilament* FilopodiaCSystem<NDIM>::initializeCFilament(int length)
    {
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
        
        ///maxlen, for now
        int maxUnits = L / monomer_size;
        int numUnits = length / monomer_size;
        
        ///initialize each filament
        CFilament* f = new CFilament();
        Compartment* cNext = cStart;
        
        int numSubFilaments = (numUnits - 1) / maxUnits + 1;
        int ci = 1;
        for(int si = 0; si < numSubFilaments; si++) {
            int setLength; ///length to intialize subfilament to
            
            if (si == numSubFilaments - 1) 
                setLength = numUnits - (maxUnits * (numSubFilaments - 1));
            else
                setLength = maxUnits;
            
            CSubFilament* currentSubFilament =
                CSystem<NDIM>::_initializer->createCSubFilament(f, cNext, {"Actin"}, setLength);
            currentSubFilament->setLength(setLength);
            
            
            if(si * L >= cStart->sides()[0] * ci) {
                cNext = cNext->neighbours().back();
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
        CSubFilament* s1 = f->getFrontCSubFilament();
        Compartment* cCurrent = s1->compartment();
        Compartment* cNext;
        
        float compartmentLength = CSystem<NDIM>::_grid->getProtoCompartment().sides()[0];
        
        if(f->length() * monomer_size >= f->numCompartments() * compartmentLength) {
            
            if(cCurrent->neighbours().size() == 1 && cCurrent != CSystem<NDIM>::_grid->getCompartment(0)) {
                std::cout << "A Filament has reached the end of the grid. Exiting." << std::endl;
                std::cout << std::endl;
                CSystem<NDIM>::printFilaments();
                exit(EXIT_SUCCESS);
            }
            
            cNext = cCurrent->neighbours().back();
            f->increaseNumCompartments();
        }
        else cNext = cCurrent;
        
        ///Initialize new subfilament
        CSystem<NDIM>::_initializer->createCSubFilament(f, cNext, species, 1);
        
        ///Increase length
        f->increaseLength();
        
    }
    
    ///Retract the front of a CFilament
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::retractFrontOfCFilament(CFilament *f)
    {
        ///Decrease length
        f->decreaseLength();
        
        ///remove front sub filament
        f->removeCSubFilament(f->getFrontCSubFilament());
        
        ///if num sub filaments has dropped to zero, remove filament from system
        if(f->numCSubFilaments() == 0) {
            auto fUnique = static_cast<std::unique_ptr<CFilament>>(f);
            auto it = std::find(CSystem<NDIM>::_filaments.begin(),
                                CSystem<NDIM>::_filaments.end(), fUnique);
            CSystem<NDIM>::_filaments.erase(it);
            return;
        }
        
        ///Set the new front subfilament to have a front
        static_cast<CMonomerBasic*>(f->getFrontCSubFilament()->
                                    frontCMonomer())->getFront()->setN(1);
    }
  
    
    ///perform one step of retrograde flow
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::retrogradeFlow()
    {
//        ///loop through all filaments, push species back
//        for(auto &fUnique : CSystem<NDIM>::_filaments) {
//            CFilament* f = fUnique.get();
//            
//            ///set front monomer and
//            static_cast<CMonomerBasic*>(f->getFrontCSubFilament()->
//                                    frontCMonomer())->getFront()->setN(1);
//        
//        }
        
    }
    
    //Specializations
    template std::vector<ReactionBase*>
        SimpleInitializer<1>::findPolymerizationReactions(CFilament* f);
    template std::vector<ReactionBase*>
        SimpleInitializer<2>::findPolymerizationReactions(CFilament* f);
    template std::vector<ReactionBase*>
        SimpleInitializer<3>::findPolymerizationReactions(CFilament* f);
    
    template void SimpleInitializer<1>::initializeProtoCompartment(CompartmentSpatial<1>& Cproto);
    template void SimpleInitializer<2>::initializeProtoCompartment(CompartmentSpatial<2>& Cproto);
    template void SimpleInitializer<3>::initializeProtoCompartment(CompartmentSpatial<3>& Cproto);
    
    template void SimpleInitializer<1>::update(CFilament* f, ReactionBase* r);
    template void SimpleInitializer<2>::update(CFilament* f, ReactionBase* r);
    template void SimpleInitializer<3>::update(CFilament* f, ReactionBase* r);
    
    template CSubFilament*
    SimpleInitializer<1>::createCSubFilament(CFilament* parentFilament,
                                             Compartment* c,
                                             std::vector<std::string> species,
                                             int length);
    template CSubFilament*
    SimpleInitializer<2>::createCSubFilament(CFilament* parentFilament,
                                             Compartment* c,
                                             std::vector<std::string> species,
                                             int length);
    template CSubFilament*
    SimpleInitializer<3>::createCSubFilament(CFilament* parentFilament,
                                             Compartment* c,
                                             std::vector<std::string> species,
                                             int length);
    
    //Specializations
    template CFilament* FilopodiaCSystem<1>::initializeCFilament(int length);
    template CFilament* FilopodiaCSystem<2>::initializeCFilament(int length);
    template CFilament* FilopodiaCSystem<3>::initializeCFilament(int length);
    
    template void
        FilopodiaCSystem<1>::extendFrontOfCFilament(CFilament *f,
                                                    std::vector<std::string> species);
    template void
        FilopodiaCSystem<2>::extendFrontOfCFilament(CFilament *f,
                                                    std::vector<std::string> species);
    template void
        FilopodiaCSystem<3>::extendFrontOfCFilament(CFilament *f,
                                                    std::vector<std::string> species);
    
    template void FilopodiaCSystem<1>::retractFrontOfCFilament(CFilament* f);
    template void FilopodiaCSystem<2>::retractFrontOfCFilament(CFilament* f);
    template void FilopodiaCSystem<3>::retractFrontOfCFilament(CFilament* f);
    
}; //end namespace chem


