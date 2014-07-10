//
//  FilamentControllerImpl.cpp
//  CytoSim
//
//  Created by James Komianos on 7/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "FilamentControllerImpl.h"

namespace chem {

    ///Connect two filaments, back to front
    ///For this impl, only add a polymerization reaction between them
    template<size_t NDIM>
    void FilopodiaInitializer<NDIM>::connect (SubFilament* s1, SubFilament* s2)
    {
        //get all species
        Compartment* c1 = s1->compartment();
        Compartment* c2 = s2->compartment();
        
        //ends
        Species* end1 = s1->backMonomer()->species("Front");
        Species* end2 = s2->frontMonomer()->species("Front");
        //Diffusing
        Species* actinDiffusing1 = c1->findSpeciesDiffusingByName("Actin");
        Species* actinDiffusing2 = c2->findSpeciesDiffusingByName("Actin");
        ///in-filament actin
        Species* actin2 = s2->frontMonomer()->species("Actin");
        
        ///Add reactions
        ReactionBase* rPoly = c1->addInternal<Reaction,2,2>({actinDiffusing1, end1, actin2, end2}, _k_on_plus);
        ReactionBase* rDepoly = c2->addInternal<Reaction,2,2>({end2, actin2, end1, actinDiffusing2}, _k_off_plus);
        
        ///Reinit chemsim
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rPoly);
        FilamentInitializer<NDIM>::_chem.addAndActivateReaction(rDepoly);
    };


    ///Initializer, based on the given simulation
    ///@param length - starting length of the filament initialized
    ///@param maxlength - length of entire reactive filament
    template <size_t NDIM>
    SubFilament* FilopodiaInitializer<NDIM>::createSubFilament(Filament* parentFilament,
                                                         int length,
                                                         int maxlength,
                                                         Compartment* c)
    {
        SubFilament* subf = new SubFilament(c);
        parentFilament->addSubFilament(subf);
        
        ///Add monomers
        for (int index = 0; index < maxlength; index++)
        {
            ///Add species we want
            SpeciesFilament* actin = c->addSpeciesFilament("Actin", (index < length) ? 1 : 0, 1);
            SpeciesFilament* capping = c->addSpeciesFilament("Capping", 0, 1);
            SpeciesFilament* formin = c->addSpeciesFilament("X-Formin", 0, 1);
            
            ///Front and back
            SpeciesFilament* back = c->addSpeciesFilament("Back", 0, 1);
            SpeciesFilament* front = c->addSpeciesFilament("Front", 0, 1);
            
            if(index == 0 && (parentFilament->numSubFilaments() == 1)) back->setN(1);
            if(index == length - 1) front->setN(1);
            
            ///add to subfilament
            subf->addMonomer(new Monomer({actin, capping, formin, back, front}, c));
        }
        
        //Look up diffusing species in this compartment
        Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
        
        ///Add polymerization reactions
        for (int index = 0; index < maxlength; index++) {
            
            ReactionBase *r;
            
            if (index != maxlength - 1) {
                r = c->addInternal<Reaction,2,2>({subf->getMonomerSpecies(index, "Front"), actinDiffusing,
                    subf->getMonomerSpecies(index + 1, "Actin"), subf->getMonomerSpecies(index + 1, "Front")}, _k_on_plus);
            }
            ///extension callback
            else {
                r = c->addInternal<Reaction,2,0>({subf->getMonomerSpecies(index, "Front"), actinDiffusing}, _k_on_plus);
                boost::signals2::shared_connection_block
                rcb(r->connect(FilamentExtensionCallback<NDIM>(FilamentInitializer<NDIM>::_controller,
                                                               parentFilament), false));
            }
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
        }
        
        /// add depolymerization reactions (+)
        for (int index = maxlength - 1; index > 0; index--) {
            
            ReactionBase* r = c->addInternal<Reaction,2,2>({subf->getMonomerSpecies(index, "Front"),
                subf->getMonomerSpecies(index, "Actin"),
                actinDiffusing, subf->getMonomerSpecies(index - 1, "Front")}, _k_off_plus);
            
            FilamentInitializer<NDIM>::_chem.addAndActivateReaction(r);
        }
        
        return subf;
    }
    
    //specializations
    template void FilopodiaInitializer<1>::connect (SubFilament* s1, SubFilament* s2);
    template void FilopodiaInitializer<2>::connect (SubFilament* s1, SubFilament* s2);
    template void FilopodiaInitializer<3>::connect (SubFilament* s1, SubFilament* s2);
    
    template SubFilament* FilopodiaInitializer<1>::createSubFilament(Filament* parentFilament,
                                                                     int length,
                                                                     int maxlength,
                                                                     Compartment* c);
    template SubFilament* FilopodiaInitializer<2>::createSubFilament(Filament* parentFilament,
                                                                     int length,
                                                                     int maxlength,
                                                                     Compartment* c);
    template SubFilament* FilopodiaInitializer<3>::createSubFilament(Filament* parentFilament,
                                                                     int length,
                                                                     int maxlength,
                                                                     Compartment* c);
    

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
            FilamentController<NDIM>::_initializer->createSubFilament(f, length, maxLength, c_start);
            
            filaments->emplace(f);
        }
        
        return filaments;
    }

    ///Extend the front of a filament
    template <size_t NDIM>
    void FilamentControllerBasic<NDIM>::extendFrontOfFilament(Filament *f)
    {
        ///Find next compartment (1D for now)
        Compartment* cCurrent = f->getFrontSubFilament()->compartment();
        Compartment* cNext = cCurrent->neighbours().back();
        
        SubFilament* s1 = f->getFrontSubFilament();
        
        ///maxlen, for now
        int maxLength = int(static_cast<CompartmentSpatial<NDIM>*>(cCurrent)->sides()[0] / monomer_size);
        
        ///Initialize new subfilament and connect it
        SubFilament* s2 = FilamentController<NDIM>::_initializer->createSubFilament(f, 1, maxLength, cNext);
        FilamentController<NDIM>::_initializer->connect(s1,s2);
        
    }

    
    //Specializations
    
    template std::unordered_set<std::unique_ptr<Filament>>*
        FilamentControllerBasic<1>::initialize(int numFilaments, int length);
    template std::unordered_set<std::unique_ptr<Filament>>*
        FilamentControllerBasic<2>::initialize(int numFilaments, int length);
    template std::unordered_set<std::unique_ptr<Filament>>*
        FilamentControllerBasic<3>::initialize(int numFilaments, int length);
    
    template void FilamentControllerBasic<1>::extendFrontOfFilament(Filament *f);
    template void FilamentControllerBasic<2>::extendFrontOfFilament(Filament *f);
    template void FilamentControllerBasic<3>::extendFrontOfFilament(Filament *f);
    
}; //end namespace chem


