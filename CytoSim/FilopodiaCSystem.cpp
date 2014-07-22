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

    ///REACTION CALLBACKS
    
    ///Extension callback
    template<size_t NDIM>
    struct CFilamentExtensionCallback {
        
        //members
        CSystem<NDIM>* _csystem;
        CFilament* _filament;
        std::vector<std::string>* _species;
        
        ///Constructor, sets members
        CFilamentExtensionCallback(CSystem<NDIM>* csystem,
                                   CFilament* filament,
                                   std::vector<std::string>* species) :
        _csystem(csystem), _filament(filament), _species(species) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            _csystem->extendFrontOfCFilament(_filament, _species);
            _csystem->update(_filament, r);
        }
    };
    
    ///General polymerization callback
    template<size_t NDIM>
    struct CFilamentPolyCallback {
        
        //members
        CSystem<NDIM> *_csystem;
        CFilament* _filament;
        
        CFilamentPolyCallback(CSystem<NDIM>* csystem,
                              CFilament* filament) :
        _csystem(csystem), _filament(filament) {};
        
        //Callback
        void operator() (ReactionBase *r){
            _filament->increaseLength();
            _csystem->update(_filament, r);
        }
        
    };
    
    ///General depolymerization callback
    template<size_t NDIM>
    struct CFilamentDepolyCallback {
        
        //members
        CSystem<NDIM> *_csystem;
        CFilament* _filament;
        
        CFilamentDepolyCallback(CSystem<NDIM>* csystem,
                                CFilament* filament) :
        _csystem(csystem), _filament(filament) {};
        
        //Callback
        void operator() (ReactionBase *r){
            _filament->decreaseLength();
            _csystem->update(_filament, r);
        }
        
    };
    
    
    ///Find the current polymerization reactions associated with this CFilament
    template<size_t NDIM>
    std::vector<ReactionBase*>*
        SimpleInitializer<NDIM>::findPolymerizationReactions(CFilament* f)
    {
        auto polyReactions = new std::vector<ReactionBase*>();
        ///go to the front subfilament, front monomer
        CMonomer* frontCMonomer = f->getFrontCSubFilament()->frontCMonomer();
        
        ///find all reactions associated with this active end polymerization
        Species* end = frontCMonomer->getActiveEndSpecies();
        auto reactions = end->getRSpecies().ReactantReactions();
        
        for(auto it = reactions.begin(); it != reactions.end(); it++)
            if((*it)->isPolymerizationReaction())
                polyReactions->push_back(*it);
  
        return polyReactions;
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

    ///Initializer, based on the given simulation
    ///@param length - starting length of the filament initialized
    ///@param maxlength - length of entire reactive filament
    template <size_t NDIM>
    CSubFilament* SimpleInitializer<NDIM>::createCSubFilament(CFilament* parentFilament,
                                                               Compartment* c,
                                                               std::vector<std::string>* species,
                                                               int length,
                                                               int maxlength)
    {
        CSubFilament* subf = new CSubFilament(c);
        parentFilament->addCSubFilament(subf);
        subf->setLength(length);
        
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
                if (std::find(species->begin(), species->end(), "X-Formin") != species->end())
                    formin->setN(1);
                else front->setN(1);
                
            }

            if(index == 0 && (parentFilament->numCSubFilaments() == 1)) back->setN(1);
            
            ///add to subfilament
            subf->addCMonomer(new CMonomerBasic({actin, front, back, formin, capping}, c));
            subf->addCBound(new CBoundBasic({empty, myosin, ma}, c));
        }
        
        ///Callbacks needed
        CFilamentDepolyCallback<NDIM> depolyCallback =
        CFilamentDepolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,parentFilament);
        
        CFilamentPolyCallback<NDIM> polyCallback =
        CFilamentPolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem, parentFilament);
        
        CFilamentExtensionCallback<NDIM> extensionCallback =
        CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                         parentFilament,
                                         new std::vector<std::string>{"Actin"});
        
        CFilamentExtensionCallback<NDIM> extensionForminCallback =
        CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                         parentFilament,
                                         new std::vector<std::string>{"Actin", "Formin"});
        
        
        //Look up diffusing species in this compartment
        Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
        Species* cappingDiffusing = c->findSpeciesDiffusingByName("Capping");
        Species* forminDiffusing = c->findSpeciesDiffusingByName("X-Formin");
        //Species* myosinDiffusing = c->findSpeciesDiffusingByName("Myosin");
        
        
        ///Loop through all spots in subfilament
        for (int index = 0; index < maxlength; index++) {
        
            ///Monomer and bounds at current index
            //CBoundBasic *b1 = static_cast<CBoundBasic*>(subf->getCBound(index));
            CBoundBasic *b2 = static_cast<CBoundBasic*>(subf->getCBound(index+1));
            
            CMonomerBasic *m1 = static_cast<CMonomerBasic*>(subf->getCMonomer(index));
            CMonomerBasic *m2 = static_cast<CMonomerBasic*>(subf->getCMonomer(index+1));
            
            ReactionBase *rPoly, *rDepoly;
            if (index < maxlength - 1) {
                
                ///Add basic polymerization reactions
                rPoly = c->addInternal<Reaction,2,3>({m1->getFront(), actinDiffusing, m2->getActin(),
                                                  m2->getFront(), b2->getEmpty()}, _k_on_plus);
                boost::signals2::shared_connection_block rcb1(rPoly->connect(polyCallback,false));
                
                /// add basic depolymerization reactions (+)
                rDepoly = c->addInternal<Reaction,3,2>({m2->getFront(), m2->getActin(), b2->getEmpty(),
                    m1->getFront(), actinDiffusing}, _k_off_plus);
                
                boost::signals2::shared_connection_block rcb2(rDepoly->connect(depolyCallback,false));
                
                CFilamentInitializer<NDIM>::_chem.addReaction(rDepoly);
            }
            ///extension callback
            else {
                rPoly = c->addInternal<Reaction,2,0>({m1->getFront(), actinDiffusing},_k_on_plus);
                
                boost::signals2::shared_connection_block rcb(rPoly->connect(extensionCallback,false));
            }
            
            rPoly->setAsPolymerizationReaction();
            CFilamentInitializer<NDIM>::_chem.addReaction(rPoly);
        
            ///add capping polymerization and depolymerization reactions (+)
            rPoly = c->addInternal<Reaction,2,1>({cappingDiffusing, m1->getFront(),
                                              m1->getCapping()}, _k_capping_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getCapping(), m1->getFront(),
                                              cappingDiffusing}, _k_capping_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            CFilamentInitializer<NDIM>::_chem.addReaction(rPoly);
            CFilamentInitializer<NDIM>::_chem.addReaction(rDepoly);
        
            ///add formin polymerization and depolymerization reactions (+)
            
            rPoly = c->addInternal<Reaction,2,1>({forminDiffusing, m1->getFront(),
                                                m1->getFormin()}, _k_formin_on_plus);
            
            rDepoly = c->addInternal<Reaction,1,2>({m1->getFormin(), m1->getFront(),
                                                forminDiffusing}, _k_formin_off_plus);
            
            rPoly->setAsPolymerizationReaction();
            CFilamentInitializer<NDIM>::_chem.addReaction(rPoly);
            CFilamentInitializer<NDIM>::_chem.addReaction(rDepoly);
            
            ///add accelerated polymerization of actin with anti-cappped end
            if (index < maxlength - 1) {
                rPoly = c->addInternal<Reaction,2,3>({actinDiffusing, m1->getFormin(), m2->getActin(),
                                                  m2->getFormin(), b2->getEmpty()}, _k_accel_on_plus);
                
                boost::signals2::shared_connection_block rcb(rPoly->connect(polyCallback, false));
            }
            ///extension callback
            else {
                rPoly = c->addInternal<Reaction,2,0>({m1->getFormin(), actinDiffusing},
                                                  _k_accel_on_plus);
            
                boost::signals2::shared_connection_block rcb(rPoly->connect(extensionForminCallback,false));
            }
            
            rPoly->setAsPolymerizationReaction();
            CFilamentInitializer<NDIM>::_chem.addReaction(rPoly);
        
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
        ///clean up and return
        delete species;
        return subf;
    }
    
    ///Connect two CFilaments, back to front
    ///For this impl, only add a polymerization reaction between them
    template<size_t NDIM>
    void SimpleInitializer<NDIM>::connect (CSubFilament* s1, CSubFilament* s2)
    {
        //get all species
        Compartment* c1 = s1->compartment();
        Compartment* c2 = s2->compartment();
        
        CFilament* parentFilament = static_cast<CFilament*>(s1->getParent());
        
        ///Monomer and bounds at current index
        //CBoundBasic *b1 = static_cast<CBoundBasic*>(s1->frontCBound());
        CBoundBasic *b2 = static_cast<CBoundBasic*>(s2->backCBound());
        
        CMonomerBasic *m1 = static_cast<CMonomerBasic*>(s1->backCMonomer());
        CMonomerBasic *m2 = static_cast<CMonomerBasic*>(s2->frontCMonomer());
        
        Species* actinDiffusing1 = c1->findSpeciesDiffusingByName("Actin");
        Species* actinDiffusing2 = c2->findSpeciesDiffusingByName("Actin");
        
        
        CFilamentPolyCallback<NDIM> polyCallback =
            CFilamentPolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                                                parentFilament);
        CFilamentDepolyCallback<NDIM> depolyCallback =
            CFilamentDepolyCallback<NDIM>(CFilamentInitializer<NDIM>::_csystem,
                                                                parentFilament);
        
        ///Add reactions
        ReactionBase* rAccelPoly =
            c1->addInternal<Reaction,2,3>({actinDiffusing1, m1->getFormin(), m2->getActin(),
                                           m2->getFormin(), b2->getEmpty()}, _k_accel_on_plus);
        boost::signals2::shared_connection_block
        rcb1(rAccelPoly->connect(polyCallback, false));
        
        ReactionBase* rPoly =
            c1->addInternal<Reaction,2,3>({actinDiffusing1, m1->getFront(), m2->getActin(),
                                           m2->getFront(), b2->getEmpty()}, _k_on_plus);
        boost::signals2::shared_connection_block
        rcb2(rPoly->connect(polyCallback, false));
        
        ReactionBase* rDepoly =
            c2->addInternal<Reaction,3,2>({m2->getFront(), m2->getActin(), b2->getEmpty(),
                                           m1->getFront(), actinDiffusing2}, _k_off_plus);
        boost::signals2::shared_connection_block
        rcb3(rDepoly->connect(depolyCallback, false));

        CFilamentInitializer<NDIM>::_chem.addReaction(rPoly);
        CFilamentInitializer<NDIM>::_chem.addReaction(rDepoly);
        CFilamentInitializer<NDIM>::_chem.addReaction(rAccelPoly);

//        ReactionBase* rMForwardStep =
//            c1->addInternal<Reaction,2,2>({b1->getMyosin(), b2->getEmpty(),
//                                           b2->getMyosin(), b1->getEmpty()}, _k_forward_step);
//        
//        ReactionBase* rMBackwardStep =
//            c2->addInternal<Reaction,2,2>({b2->getMyosin(), b1->getEmpty(),
//                                           b1->getMyosin(), b2->getEmpty()}, _k_backward_step);
//        
//        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMForwardStep);
//        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMBackwardStep);
//        
//        ReactionBase* rMA3ForwardStep =
//            c1->addInternal<Reaction,2,2>({b1->getMyosinActin(), b2->getEmpty(),
//                                           b2->getMyosinActin(), b1->getMyosinActin()},
//                                           _k_forward_step);
//        
//        ReactionBase* rMA3BackwardStep =
//            c2->addInternal<Reaction,2,2>({b2->getMyosinActin(), b1->getEmpty(),
//                                           b1->getMyosinActin(), b2->getEmpty()},
//                                           _k_backward_step);
//        
//        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3ForwardStep);
//        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3BackwardStep);
    };
    
    //specializations
    template std::vector<ReactionBase*>*
        SimpleInitializer<1>::findPolymerizationReactions(CFilament* f);
    template std::vector<ReactionBase*>*
        SimpleInitializer<2>::findPolymerizationReactions(CFilament* f);
    template std::vector<ReactionBase*>*
        SimpleInitializer<3>::findPolymerizationReactions(CFilament* f);
    
    template void SimpleInitializer<1>::update(CFilament* f, ReactionBase* r);
    template void SimpleInitializer<2>::update(CFilament* f, ReactionBase* r);
    template void SimpleInitializer<3>::update(CFilament* f, ReactionBase* r);
    
    template void SimpleInitializer<1>::connect (CSubFilament* s1, CSubFilament* s2);
    template void SimpleInitializer<2>::connect (CSubFilament* s1, CSubFilament* s2);
    template void SimpleInitializer<3>::connect (CSubFilament* s1, CSubFilament* s2);
    
    template CSubFilament*
        SimpleInitializer<1>::createCSubFilament(CFilament* parentFilament,
                                                    Compartment* c,
                                                    std::vector<std::string>* species,
                                                    int length,
                                                    int maxlength);
    template CSubFilament*
        SimpleInitializer<2>::createCSubFilament(CFilament* parentFilament,
                                                    Compartment* c,
                                                    std::vector<std::string>* species,
                                                    int length,
                                                    int maxlength);
    template CSubFilament*
        SimpleInitializer<3>::createCSubFilament(CFilament* parentFilament,
                                                    Compartment* c,
                                                    std::vector<std::string>* species,
                                                    int length,
                                                    int maxlength);
    

    //Initialize a number of filaments
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::initialize(int numFilaments, int length)
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
        int maxLength = int(cStart->sides()[0] / monomer_size);
        
        ///initialize each filament
        for(int fi = 0; fi < numFilaments; fi++) {
            
            CFilament* f = new CFilament();
            Compartment* cNext = cStart;
            CSubFilament* lastSubFilament = nullptr;
            
            int numSubFilaments = length / (maxLength + 1) + 1;
            
            for(int si = 0; si < numSubFilaments; si++) {
                int setLength; ///length to intialize subfilament to
                
                if (si == numSubFilaments - 1)
                    setLength = length - (maxLength * (numSubFilaments - 1));
                else
                    setLength = maxLength;
                
                CSubFilament* currentSubFilament =
                    CSystem<NDIM>::_initializer->createCSubFilament(f, cNext,
                             new std::vector<std::string>{"Actin"}, setLength, maxLength);
                
                if(lastSubFilament != nullptr)
                    CSystem<NDIM>::_initializer->
                                        connect(lastSubFilament, currentSubFilament);
                    
                lastSubFilament = currentSubFilament;
                cNext = cNext->neighbours().back();
            }
            f->setLength(length);
            CSystem<NDIM>::_filaments.emplace(f);
            CSystem<NDIM>::update(f, nullptr);
            
        }
    }

    ///Extend the front of a filament
    template <size_t NDIM>
    void FilopodiaCSystem<NDIM>::extendFrontOfCFilament(CFilament *f,
                                                           std::vector<std::string>* species)
    {
        ///Find next compartment (1D for now)
        Compartment* cCurrent = f->getFrontCSubFilament()->compartment();
        Compartment* cNext = cCurrent->neighbours().back();
        
        CSubFilament* s1 = f->getFrontCSubFilament();
        
        ///maxlen, for now
        int sideLength = static_cast<CompartmentSpatial<NDIM>*>(cCurrent)->sides()[0];
        int maxLength = int(sideLength / monomer_size);
        
        ///Initialize new subfilament and connect it
        CSubFilament* s2 = CSystem<NDIM>::_initializer->
                            createCSubFilament(f, cNext, species, 1, maxLength);
        CSystem<NDIM>::_initializer->connect(s1,s2);
        
        ///Increase length
        f->increaseLength();
        
    }

    
    //Specializations
    
    template void FilopodiaCSystem<1>::initialize(int numFilaments, int length);
    template void FilopodiaCSystem<2>::initialize(int numFilaments, int length);
    template void FilopodiaCSystem<3>::initialize(int numFilaments, int length);
    
    template void
        FilopodiaCSystem<1>::extendFrontOfCFilament(CFilament *f,
                                                       std::vector<std::string>* species);
    template void
        FilopodiaCSystem<2>::extendFrontOfCFilament(CFilament *f,
                                                       std::vector<std::string>* species);
    template void
        FilopodiaCSystem<3>::extendFrontOfCFilament(CFilament *f,
                                                       std::vector<std::string>* species);
    
}; //end namespace chem


