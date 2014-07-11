//
//  FilamentControllerImpl.cpp
//  CytoSim
//
//  Created by James Komianos on 7/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "FilamentControllerImpl.h"

namespace chem {

    ///Initializer, based on the given simulation
    ///@param length - starting length of the filament initialized
    ///@param maxlength - length of entire reactive filament
    template <size_t NDIM>
    SubFilament* FilopodiaInitializer<NDIM>::createSubFilament(Filament* parentFilament,
                                                               Compartment* c,
                                                               std::vector<std::string>* species,
                                                               int length,
                                                               int maxlength)
    {
        SubFilament* subf = new SubFilament(c);
        parentFilament->addSubFilament(subf);
        
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

            if(index == 0 && (parentFilament->numSubFilaments() == 1)) back->setN(1);
            
            ///add to subfilament
            subf->addMonomer(new Monomer({actin, capping, formin, back, front}, c));
            subf->addBound(new Bound({myosin, ma, empty}, c));
        }
        
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
            }
            ///extension callback
            else {
                r = c->addInternal<Reaction,2,0>({subf->getMonomerSpecies(index, "Front"),
                                                  actinDiffusing},
                                                  _k_on_plus);
                
                boost::signals2::shared_connection_block
                rcb(r->connect(FilamentExtensionCallback<NDIM>(FilamentInitializer<NDIM>::_controller,
                                                               parentFilament,
                                                               new std::vector<std::string>{"Actin"}),
                                                               false));
            }
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
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
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
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
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
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
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
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
            }
            ///extension callback
            else {
                r = c->addInternal<Reaction,2,0>({subf->getMonomerSpecies(index, "X-Formin"),
                                                  actinDiffusing},
                                                  _k_accel_on_plus);
            
                boost::signals2::shared_connection_block
                rcb(r->connect(FilamentExtensionCallback<NDIM>(FilamentInitializer<NDIM>::_controller,
                                                               parentFilament,
                                                               new std::vector<std::string>{"Actin","X-Formin"}),
                                                               false));
            }
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
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
            
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rBinding);
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rUnbinding);
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rLoading);
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rUnloading);

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
                
                FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMForwardStep);
            }
            if (index != 0) {
                rMBackwardStep =
                    c->addInternal<Reaction,2,2>({subf->getBoundSpecies(index, "Myosin"),
                                                  subf->getBoundSpecies(index-1, "Empty"),
                                                  subf->getBoundSpecies(index-1, "Myosin"),
                                                  subf->getBoundSpecies(index, "Empty")},
                                                  _k_forward_step);
                
                FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMBackwardStep);
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
                
                FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3ForwardStep);
            }
            if (index != 0) {
                rMA3BackwardStep =
                    c->addInternal<Reaction,2,2>({subf->getBoundSpecies(index, "A-MyosinActin"),
                                                  subf->getBoundSpecies(index-1, "Empty"),
                                                  subf->getBoundSpecies(index-1, "A-MyosinActin"),
                                                  subf->getBoundSpecies(index, "Empty")},
                                                  _k_forward_step);
                
                FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3BackwardStep);
            }
            
            ReactionBase *rMA3Unbinding =
                c->addInternal<Reaction,1,3>({subf->getBoundSpecies(index, "A-MyosinActin"),
                                              subf->getBoundSpecies(index, "Myosin"),
                                              subf->getBoundSpecies(index, "Empty"),
                                              actinDiffusing},
                                              _k_unbinding);
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3Unbinding);
        }
        
        ///clean up and return
        delete species;
        return subf;
    }
    
    ///Connect two filaments, back to front
    ///For this impl, only add a polymerization reaction between them
    template<size_t NDIM>
    void FilopodiaInitializer<NDIM>::connect (SubFilament* s1, SubFilament* s2)
    {
        //get all species
        Compartment* c1 = s1->compartment();
        Compartment* c2 = s2->compartment();
        
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
        ReactionBase* rPoly =
            c1->addInternal<Reaction,2,3>({actinDiffusing1, front1,
                                           actin2, front2, empty2},
                                           _k_on_plus);
        ReactionBase* rDepoly =
            c2->addInternal<Reaction,3,2>({front2, actin2, empty2,
                                           front1, actinDiffusing2},
                                           _k_off_plus);
        
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
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rAccelPoly);
        
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMForwardStep);
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMBackwardStep);

        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3ForwardStep);
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rMA3BackwardStep);
    };

    
    
    //specializations
    template void FilopodiaInitializer<1>::connect (SubFilament* s1, SubFilament* s2);
    template void FilopodiaInitializer<2>::connect (SubFilament* s1, SubFilament* s2);
    template void FilopodiaInitializer<3>::connect (SubFilament* s1, SubFilament* s2);
    
    template SubFilament* FilopodiaInitializer<1>::createSubFilament(Filament* parentFilament,
                                                                     Compartment* c,
                                                                     std::vector<std::string>* species,
                                                                     int length,
                                                                     int maxlength);
    template SubFilament* FilopodiaInitializer<2>::createSubFilament(Filament* parentFilament,
                                                                     Compartment* c,
                                                                     std::vector<std::string>* species,
                                                                     int length,
                                                                     int maxlength);
    template SubFilament* FilopodiaInitializer<3>::createSubFilament(Filament* parentFilament,
                                                                     Compartment* c,
                                                                     std::vector<std::string>* species,
                                                                     int length,
                                                                     int maxlength);
    

    //Initialize a number of filaments
    template <size_t NDIM>
    std::unordered_set<std::unique_ptr<Filament>>* FilamentControllerBasic<NDIM>::initialize(int numFilaments, int length)
    {
        CompartmentSpatial<NDIM>* c_start;
        ///Starting compartment for 1D, all filaments start in compartment 0
        if (NDIM == 1) {
            c_start = FilamentController<NDIM>::_grid->getCompartment(0);
        }
        else {
            std::cout << "Multiple dimensional implementation not done. Exiting." << std::endl;
            exit(EXIT_FAILURE);
        }
        ///maxlen, for now
        int maxLength = int(c_start->sides()[0] / monomer_size);
        
        ///set of filaments
        auto filaments = new std::unordered_set<std::unique_ptr<Filament>>() ;
        
        ///initialize each filament
        for(int index = 0; index < numFilaments; index++) {
            
            Filament* f = new Filament();
            FilamentController<NDIM>::_initializer->createSubFilament(f, c_start,
                            new std::vector<std::string>{"Actin"}, length, maxLength);
            
            filaments->emplace(f);
        }
        
        return filaments;
    }

    ///Extend the front of a filament
    template <size_t NDIM>
    void FilamentControllerBasic<NDIM>::extendFrontOfFilament(Filament *f, std::vector<std::string>* species)
    {
        ///Find next compartment (1D for now)
        Compartment* cCurrent = f->getFrontSubFilament()->compartment();
        Compartment* cNext = cCurrent->neighbours().back();
        
        SubFilament* s1 = f->getFrontSubFilament();
        
        ///maxlen, for now
        int maxLength = int(static_cast<CompartmentSpatial<NDIM>*>(cCurrent)->sides()[0] / monomer_size);
        
        ///Initialize new subfilament and connect it
        SubFilament* s2 = FilamentController<NDIM>::_initializer->createSubFilament(f, cNext, species, 1, maxLength);
        FilamentController<NDIM>::_initializer->connect(s1,s2);
        
    }

    
    //Specializations
    
    template std::unordered_set<std::unique_ptr<Filament>>*
        FilamentControllerBasic<1>::initialize(int numFilaments, int length);
    template std::unordered_set<std::unique_ptr<Filament>>*
        FilamentControllerBasic<2>::initialize(int numFilaments, int length);
    template std::unordered_set<std::unique_ptr<Filament>>*
        FilamentControllerBasic<3>::initialize(int numFilaments, int length);
    
    template void FilamentControllerBasic<1>::extendFrontOfFilament(Filament *f,
                                                                    std::vector<std::string>* species);
    template void FilamentControllerBasic<2>::extendFrontOfFilament(Filament *f,
                                                                    std::vector<std::string>* species);
    template void FilamentControllerBasic<3>::extendFrontOfFilament(Filament *f,
                                                                    std::vector<std::string>* species);
    
}; //end namespace chem


