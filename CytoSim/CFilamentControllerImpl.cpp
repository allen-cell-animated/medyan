//
//  CFilamentControllerImpl.cpp
//  CytoSim
//
//  Created by James Komianos on 7/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CFilamentControllerImpl.h"

namespace chem {

    ///Initializer, based on the given simulation
    ///@param length - starting length of the filament initialized
    ///@param maxlength - length of entire reactive filament
    template <size_t NDIM>
    CSubFilament* FilopodiaInitializer<NDIM>::createCSubFilament(CFilament* parentFilament,
                                                               Compartment* c,
                                                               std::vector<std::string>* species,
                                                               int length,
                                                               int maxlength)
    {
        CSubFilament* subf = new CSubFilament(c);
        parentFilament->addCSubFilament(subf);
        
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
            
            
            if(index < length) {
                actin->setN(1);
                empty->setN(1);
            }
            if(index == length - 1) {
                if (std::find(species->begin(), species->end(), "X-Formin") != species->end())
                    formin->setN(1);
                else front->setN(1);
                
            }

            if(index == 0 && (parentFilament->numCSubFilaments() == 1)) back->setN(1);
            
            ///add to subfilament
            subf->addMonomer(new Monomer({actin, capping, formin, back, front}, c));
            subf->addBound(new Bound({myosin, ma, empty}, c));
        }
        
        ///Callbacks needed
        CFilamentPolyCallback<NDIM> polyCallback = CFilamentPolyCallback<NDIM>(parentFilament);
        
        CFilamentDepolyCallback<NDIM> depolyCallback = CFilamentDepolyCallback<NDIM>(parentFilament);
        
        CFilamentExtensionCallback<NDIM> extensionCallback =
                CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_controller,
                                                 parentFilament,
                                                 new std::vector<std::string>{"Actin"});
        CFilamentExtensionCallback<NDIM> extensionForminCallback =
            CFilamentExtensionCallback<NDIM>(CFilamentInitializer<NDIM>::_controller,
                                             parentFilament,
                                             new std::vector<std::string>{"Actin", "Formin"});
        
        
        //Look up diffusing species in this compartment
        Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
        
        ///Add basic polymerization reactions
        for (int index = 0; index < maxlength; index++) {
            
            ReactionBase *r;
            
            if (index != maxlength - 1) {
                r = c->addInternal<Reaction,2,3>({subf->getMonomerSpecies(index, "Front"),
                                                  actinDiffusing,
                                                  subf->getMonomerSpecies(index + 1, "Actin"),
                                                  subf->getMonomerSpecies(index + 1, "Front"),
                                                  subf->getBoundSpecies(index + 1, "Empty")},
                                                  _k_on_plus);
                boost::signals2::shared_connection_block rcb(r->connect(polyCallback,false));

            }
            ///extension callback
            else {
                r = c->addInternal<Reaction,2,0>({subf->getMonomerSpecies(index, "Front"),
                                                  actinDiffusing},
                                                  _k_on_plus);
                
                boost::signals2::shared_connection_block rcb(r->connect(extensionCallback,false));
            }
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
        }
        
        /// add basic depolymerization reactions (+)
        for (int index = maxlength - 1; index > 0; index--) {
            
            ReactionBase* r =
                c->addInternal<Reaction,3,2>({subf->getMonomerSpecies(index, "Front"),
                                              subf->getMonomerSpecies(index, "Actin"),
                                              subf->getBoundSpecies(index, "Empty"),
                                              subf->getMonomerSpecies(index - 1, "Front"),
                                              actinDiffusing},
                                              _k_off_plus);
            
            boost::signals2::shared_connection_block rcb(r->connect(depolyCallback,false));
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
        }
        
        Species* cappingDiffusing = c->findSpeciesDiffusingByName("Capping");
        
        ///add capping polymerization and depolymerization reactions (+)
        for (int index = 0; index < maxlength; index++) {

            ReactionBase *rPoly =
                c->addInternal<Reaction,2,1>({cappingDiffusing,
                                              subf->getMonomerSpecies(index, "Front"),
                                              subf->getMonomerSpecies(index, "Capping")},
                                              _k_capping_on_plus);
            
            ReactionBase *rDepoly =
                c->addInternal<Reaction,1,2>({subf->getMonomerSpecies(index, "Capping"),
                                              subf->getMonomerSpecies(index, "Front"),
                                              cappingDiffusing},
                                              _k_capping_off_plus);
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
        }
        
        Species* forminDiffusing = c->findSpeciesDiffusingByName("X-Formin");
        
        ///add formin polymerization and depolymerization reactions (+)
        for (int index = 0; index < maxlength; index++) {
            
            ReactionBase *rPoly =
                c->addInternal<Reaction,2,1>({forminDiffusing,
                                              subf->getMonomerSpecies(index, "Front"),
                                              subf->getMonomerSpecies(index, "X-Formin")},
                                              _k_formin_on_plus);
            
            ReactionBase *rDepoly =
                c->addInternal<Reaction,1,2>({subf->getMonomerSpecies(index, "X-Formin"),
                                              subf->getMonomerSpecies(index, "Front"),
                                              forminDiffusing},
                                               _k_formin_off_plus);
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
        }
        
        ///add accelerated polymerization of actin with anti-cappped end
        for (int index = 0; index < maxlength; index++) {
            
            ReactionBase *r;
            
            if (index != maxlength - 1) {
                r = c->addInternal<Reaction,2,3>({actinDiffusing,
                                                  subf->getMonomerSpecies(index, "X-Formin"),
                                                  subf->getMonomerSpecies(index + 1, "Actin"),
                                                  subf->getMonomerSpecies(index + 1, "X-Formin"),
                                                  subf->getBoundSpecies(index + 1, "Empty")},
                                                  _k_accel_on_plus);
                boost::signals2::shared_connection_block rcb(r->connect(polyCallback, false));
            }
            ///extension callback
            else {
                r = c->addInternal<Reaction,2,0>({subf->getMonomerSpecies(index, "X-Formin"),
                                                  actinDiffusing},
                                                  _k_accel_on_plus);
            
                boost::signals2::shared_connection_block rcb(r->connect(extensionForminCallback,false));
            }
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
        }
        
        
        ///add motor binding and unbinding, loading and unloading
        Species* myosinDiffusing = c->findSpeciesDiffusingByName("Myosin");
        
        for (int index = 0; index < maxlength; index++) {
            
            ReactionBase *rBinding =
                c->addInternal<Reaction,2,1>({myosinDiffusing,
                                              subf->getBoundSpecies(index, "Empty"),
                                              subf->getBoundSpecies(index, "Myosin")},
                                              _k_binding);
            
            ReactionBase *rUnbinding =
                c->addInternal<Reaction,1,2>({subf->getBoundSpecies(index, "Myosin"),
                                              subf->getBoundSpecies(index, "Empty"),
                                              myosinDiffusing},
                                              _k_unbinding);
            
            ReactionBase *rLoading =
                c->addInternal<Reaction,2,1>({actinDiffusing,
                                              subf->getBoundSpecies(index, "Myosin"),
                                              subf->getBoundSpecies(index, "A-MyosinActin")},
                                              _k_load);
            
            ReactionBase *rUnloading =
                c->addInternal<Reaction,1,2>({subf->getBoundSpecies(index, "A-MyosinActin"),
                                              subf->getBoundSpecies(index, "Myosin"),
                                              actinDiffusing},
                                              _k_unload);
            
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rBinding);
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rUnbinding);
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rLoading);
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rUnloading);

            //add motor stepping
            ReactionBase *rMForwardStep;
            ReactionBase *rMBackwardStep;
            
            if (index != maxlength - 1) {
                rMForwardStep =
                    c->addInternal<Reaction,2,2>({subf->getBoundSpecies(index, "Myosin"),
                                                  subf->getBoundSpecies(index+1, "Empty"),
                                                  subf->getBoundSpecies(index+1, "Myosin"),
                                                  subf->getBoundSpecies(index, "Empty")},
                                                  _k_forward_step);
                
                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMForwardStep);
            }
            if (index != 0) {
                rMBackwardStep =
                    c->addInternal<Reaction,2,2>({subf->getBoundSpecies(index, "Myosin"),
                                                  subf->getBoundSpecies(index-1, "Empty"),
                                                  subf->getBoundSpecies(index-1, "Myosin"),
                                                  subf->getBoundSpecies(index, "Empty")},
                                                  _k_forward_step);
                
                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMBackwardStep);
            }
            
            ReactionBase *rMA3ForwardStep;
            ReactionBase *rMA3BackwardStep;
            if (index != maxlength - 1) {
                rMA3ForwardStep =
                    c->addInternal<Reaction,2,2>({subf->getBoundSpecies(index, "A-MyosinActin"),
                                                  subf->getBoundSpecies(index+1, "Empty"),
                                                  subf->getBoundSpecies(index+1, "A-MyosinActin"),
                                                  subf->getBoundSpecies(index, "Empty")},
                                                  _k_forward_step);
                
                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3ForwardStep);
            }
            if (index != 0) {
                rMA3BackwardStep =
                    c->addInternal<Reaction,2,2>({subf->getBoundSpecies(index, "A-MyosinActin"),
                                                  subf->getBoundSpecies(index-1, "Empty"),
                                                  subf->getBoundSpecies(index-1, "A-MyosinActin"),
                                                  subf->getBoundSpecies(index, "Empty")},
                                                  _k_forward_step);
                
                CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3BackwardStep);
            }
            
            ReactionBase *rMA3Unbinding =
                c->addInternal<Reaction,1,3>({subf->getBoundSpecies(index, "A-MyosinActin"),
                                              subf->getBoundSpecies(index, "Myosin"),
                                              subf->getBoundSpecies(index, "Empty"),
                                              actinDiffusing},
                                              _k_unbinding);
            
            CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3Unbinding);
        }
        
        ///clean up and return
        delete species;
        return subf;
    }
    
    ///Connect two CFilaments, back to front
    ///For this impl, only add a polymerization reaction between them
    template<size_t NDIM>
    void FilopodiaInitializer<NDIM>::connect (CSubFilament* s1, CSubFilament* s2)
    {
        //get all species
        Compartment* c1 = s1->compartment();
        Compartment* c2 = s2->compartment();
        
        CFilament* parentFilament = static_cast<CFilament*>(s1->getParent());
        
        //Ends
        Species* front1 = s1->backMonomer()->species("Front");
        Species* front2 = s2->frontMonomer()->species("Front");
        
        //Formin
        Species* formin1 = s1->backMonomer()->species("X-Formin");
        Species* formin2 = s2->frontMonomer()->species("X-Formin");
        
        //Diffusing
        Species* actinDiffusing1 = c1->findSpeciesDiffusingByName("Actin");
        Species* actinDiffusing2 = c2->findSpeciesDiffusingByName("Actin");
        
        ///in-filament actin
        Species* actin2 = s2->frontMonomer()->species("Actin");
        
        ///motors
        Species* myosin1 = s1->backBound()->species("Myosin");
        Species* myosin2 = s2->frontBound()->species("Myosin");
        
        Species* ma1 = s1->backBound()->species("A-MyosinActin");
        Species* ma2 = s2->frontBound()->species("A-MyosinActin");
        
        ///emptys
        Species* empty1 = s1->backBound()->species("Empty");
        Species* empty2 = s2->frontBound()->species("Empty");
        
        ///Add reactions 
        ReactionBase* rAccelPoly =
            c1->addInternal<Reaction,2,3>({actinDiffusing1, formin1,
                                           actin2, formin2, empty2},
                                           _k_accel_on_plus);
        boost::signals2::shared_connection_block
        rcb1(rAccelPoly->connect(CFilamentPolyCallback<NDIM>(parentFilament),
                                                   false));
        
        ReactionBase* rPoly =
            c1->addInternal<Reaction,2,3>({actinDiffusing1, front1,
                                           actin2, front2, empty2},
                                           _k_on_plus);
        boost::signals2::shared_connection_block
        rcb2(rPoly->connect(CFilamentPolyCallback<NDIM>(parentFilament),
                                                       false));
        
        ReactionBase* rDepoly =
            c2->addInternal<Reaction,3,2>({front2, actin2, empty2,
                                           front1, actinDiffusing2},
                                           _k_off_plus);
        boost::signals2::shared_connection_block
        rcb3(rDepoly->connect(CFilamentDepolyCallback<NDIM>(parentFilament),
                                                           false));
        
        
        ReactionBase* rMForwardStep =
            c1->addInternal<Reaction,2,2>({myosin1, empty2,
                                           myosin2, empty1},
                                           _k_forward_step);
        
        ReactionBase* rMBackwardStep =
            c2->addInternal<Reaction,2,2>({myosin2, empty1,
                                           myosin1, empty2},
                                           _k_backward_step);
        
        ReactionBase* rMA3ForwardStep =
            c1->addInternal<Reaction,2,2>({ma1, empty2,
                                           ma2, empty1},
                                           _k_forward_step);
        
        ReactionBase* rMA3BackwardStep =
            c2->addInternal<Reaction,2,2>({ma2, empty1,
                                           ma1, empty2},
                                           _k_backward_step);
        
        ///Reinit chemsim
        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rAccelPoly);
        
        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMForwardStep);
        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMBackwardStep);

        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3ForwardStep);
        CFilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3BackwardStep);
    };

    
    
    //specializations
    template void FilopodiaInitializer<1>::connect (CSubFilament* s1, CSubFilament* s2);
    template void FilopodiaInitializer<2>::connect (CSubFilament* s1, CSubFilament* s2);
    template void FilopodiaInitializer<3>::connect (CSubFilament* s1, CSubFilament* s2);
    
    template CSubFilament* FilopodiaInitializer<1>::createCSubFilament(CFilament* parentFilament,
                                                                     Compartment* c,
                                                                     std::vector<std::string>* species,
                                                                     int length,
                                                                     int maxlength);
    template CSubFilament* FilopodiaInitializer<2>::createCSubFilament(CFilament* parentFilament,
                                                                     Compartment* c,
                                                                     std::vector<std::string>* species,
                                                                     int length,
                                                                     int maxlength);
    template CSubFilament* FilopodiaInitializer<3>::createCSubFilament(CFilament* parentFilament,
                                                                     Compartment* c,
                                                                     std::vector<std::string>* species,
                                                                     int length,
                                                                     int maxlength);
    

    //Initialize a number of filaments
    template <size_t NDIM>
    void CFilamentControllerFilopodia<NDIM>::initialize(int numFilaments, int length)
    {
        
        ///init filament container
        CFilamentController<NDIM>::_filaments = new std::unordered_set<std::unique_ptr<CFilament>>();
        
        CompartmentSpatial<NDIM>* c_start;
        ///Starting compartment for 1D, all filaments start in compartment 0
        if (NDIM == 1) {
            c_start = CFilamentController<NDIM>::_grid->getCompartment(0);
        }
        else {
            std::cout << "Multiple dimensional implementation not done. Exiting." << std::endl;
            exit(EXIT_FAILURE);
        }
        ///maxlen, for now
        int maxLength = int(c_start->sides()[0] / monomer_size);
        
        ///initialize each filament
        for(int index = 0; index < numFilaments; index++) {
            
            CFilament* f = new CFilament();
            CFilamentController<NDIM>::_initializer->createCSubFilament(f, c_start,
                            new std::vector<std::string>{"Actin"}, length, maxLength);
            
            CFilamentController<NDIM>::_filaments->emplace(f);
            f->setLength(length);
        }
    }

    ///Extend the front of a filament
    template <size_t NDIM>
    void CFilamentControllerFilopodia<NDIM>::extendFrontOfCFilament(CFilament *f, std::vector<std::string>* species)
    {
        ///Find next compartment (1D for now)
        Compartment* cCurrent = f->getFrontCSubFilament()->compartment();
        Compartment* cNext = cCurrent->neighbours().back();
        
        CSubFilament* s1 = f->getFrontCSubFilament();
        
        ///maxlen, for now
        int maxLength = int(static_cast<CompartmentSpatial<NDIM>*>(cCurrent)->sides()[0] / monomer_size);
        
        ///Initialize new subfilament and connect it
        CSubFilament* s2 = CFilamentController<NDIM>::_initializer->
                                        createCSubFilament(f, cNext, species, 1, maxLength);
        CFilamentController<NDIM>::_initializer->connect(s1,s2);
        
        ///Increase length
        f->increaseLength();
        
    }

    
    //Specializations
    
    template void CFilamentControllerFilopodia<1>::initialize(int numFilaments, int length);
    template void CFilamentControllerFilopodia<2>::initialize(int numFilaments, int length);
    template void CFilamentControllerFilopodia<3>::initialize(int numFilaments, int length);
    
    template void CFilamentControllerFilopodia<1>::extendFrontOfCFilament(CFilament *f,
                                                                    std::vector<std::string>* species);
    template void CFilamentControllerFilopodia<2>::extendFrontOfCFilament(CFilament *f,
                                                                    std::vector<std::string>* species);
    template void CFilamentControllerFilopodia<3>::extendFrontOfCFilament(CFilament *f,
                                                                    std::vector<std::string>* species);
    
}; //end namespace chem


