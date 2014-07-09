//
//  FilamentController.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__FilamentController__
#define __CytoSim__FilamentController__

#include <iostream>
#include "Filament.h"

namespace chem {
    
    //Forward declaration
    template<size_t NDIM>
    class FilamentController;
    
    /// FilamentInitailizer class is used for initailizing a subfilament in a compartment
    /*!
     *  FilamentInitializer is an abstract class used for initailizing a subfilament, including its
     *  associated species and reactions. SubFilaments can be initialized by constructing
     *  this object and calling its initializer.
     */
    template<size_t NDIM>
    class FilamentInitializer {
        
    protected:
        FilamentController<NDIM>* _controller; ///< ptr to controller
        
    public:
        ///Constructor, destructor
        FilamentInitializer() {}
        virtual ~FilamentInitializer() {};
        
        ///Set the filament controller
        virtual void setFilamentController (FilamentController<NDIM>* controller)
        { FilamentInitializer<NDIM>::_controller = controller;}
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the SubFilament initialized
        virtual SubFilament* create(int length, int maxlength, Compartment* c) = 0;
        
        ///Connect two filaments, back to front
        virtual void connect(SubFilament* s1, SubFilament* s2) = 0;
    };
    
    
    /// FilamentController class is an abstract class for the updating of filaments
    /*!
     *  FilamentController is an abstract class that provides methods to control the updating of
     *  filaments after reactions occur. This includes setting up reactions, updating species
     *  and compartments, and creating new subfilaments.
     */
    template<size_t NDIM>
    class FilamentController {

    protected:
        CompartmentGrid<NDIM>* _grid; ///< ptr to compartment grid for updating
        FilamentInitializer<NDIM>* _initializer; ///<ptr to initializer, could be any implementation
        
    public:
        
        ///constructor and destructor
        FilamentController(CompartmentGrid<NDIM>* grid, FilamentInitializer<NDIM>* initializer)
        : _grid(grid), _initializer(initializer)
        {
            _initializer->setFilamentController(this);
        }
        
        virtual ~FilamentController() {}
        
        ///Initialize all filaments and reactions, return set of filaments
        virtual std::unordered_set<std::unique_ptr<Filament>>* initialize(int numFilaments, int length) = 0;
        
        ///Extend the front of a filament
        virtual void extendFrontOfFilament(Filament *f) = 0;
        
//        ///Retract the front of a filament
//        virtual void retractFrontOfFilament(Filament *f) = 0;
//        
//        ///Retract the back of a filament
//        virtual void retractBackOfFilament(Filament *f) = 0;
        
    };
    
    ///FILAMENT REACTION CALLBACKS
    
    template<size_t NDIM>
    struct FilamentExtensionCallback {
        
        //members
        FilamentController<NDIM>* _controller;
        Filament* _filament;
        
        
        ///Constructor, sets members
        FilamentExtensionCallback(FilamentController<NDIM>* controller, Filament* filament) :
            _controller(controller), _filament(filament) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            _controller->extendFrontOfFilament(_filament);
        }
    };
    
    /// ActinBasic is a basic implementation of the Initializer class that only involves the following:
    /// Polymerization of + end
    template<size_t NDIM>
    class ActinBasicInitializer : public FilamentInitializer<NDIM> {
        
    private:
        ///REACTION RATES
        float _k_on_plus = 0;
        float _k_off_plus = 0;
        float _k_off_minus = 0;
        
    public:
        
        ///Constructor, does nothing
        ActinBasicInitializer() {};
        
        ///Destructor, does nothing
        ~ActinBasicInitializer() {};
        
        ///Set reaction rates for this subfilament
        virtual void setReactionRates(float kOnPlus = 21.0,
                                      float kOffPlus = 1.4,
                                      float kOffMinus = 1.4)
        {
            _k_on_plus = kOnPlus;
            _k_off_plus = kOffPlus;
            _k_off_minus = kOffMinus;
        }
        
        ///Connect two filaments, back to front
        ///For this impl, only add a polymerization reaction between them
        virtual void connect (SubFilament* s1, SubFilament* s2)
        {
            //get all species
            Compartment* c = s1->compartment();

            Species* end1 = s1->backMonomer()->species("Front");
            Species* actin1 = c->findSpeciesDiffusingByName("Actin");
            Species* actin2 = s2->frontMonomer()->species("Actin");
            Species* end2 = s2->frontMonomer()->species("Front");
            
            c->addInternal<Reaction,2,2>({actin1, end1, actin2, end2}, _k_on_plus);
        };
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the filament initialized
        ///@param maxlength - length of entire reactive filament 
        virtual SubFilament* create(int length, int maxlength, Compartment* c)
        {
            SubFilament* subf = new SubFilament(c);
            
            ///Add monomers
            for (int index = 0; index < maxlength; index++)
            {
                ///Add species we want
                SpeciesFilament* actin = c->addSpeciesFilament("Actin", (index < length) ? 1 : 0);
                
                ///Front and back
                SpeciesFilament* back = c->addSpeciesFilament("Back", 0);
                SpeciesFilament* front = c->addSpeciesFilament("Front", 0);
                if(index == 0) back->setN(1);
                if(index == length - 1) front->setN(1);
                
                ///add to subfilament
                subf->addMonomer(new Monomer({actin, back, front}, c));
            }
            
            //Look up diffusing species in this compartment
            Species* actinDiffusing = c->findSpeciesDiffusingByName("Actin");
            
            ///Add polymerization reactions
            for (int index = 0; index < length; index++) {
                
                if (index != length - 1) {
                    c->addInternal<Reaction,2,2>({subf->getMonomerSpecies(index, "Front"), actinDiffusing,
                    subf->getMonomerSpecies(index + 1, "Actin"), subf->getMonomerSpecies(index + 1, "Front")}, _k_on_plus);
                }
                ///extension callback
                else {
                    ReactionBase* r =
                    c->addInternal<Reaction,2,0>({subf->getMonomerSpecies(index, "Front"), actinDiffusing}, _k_on_plus);
                    boost::signals2::shared_connection_block
                    rcb(r->connect(FilamentExtensionCallback<NDIM>(
                            FilamentInitializer<NDIM>::_controller, static_cast<Filament*>(subf->getParent())), false));
                }
            }
            return subf;
        }
        
    };
    

    /// FilamentControllerBasic is a basic implementation for updating filaments
    template <size_t NDIM>
    class FilamentControllerBasic : public FilamentController<NDIM> {
        
    public:
        
        ///Constructor, calls base class
        FilamentControllerBasic(CompartmentGrid<NDIM>* grid, FilamentInitializer<NDIM>* initializer)
        : FilamentController<NDIM>::FilamentController(grid, initializer) {};
        
        ///Default destructor, does nothing
        ~FilamentControllerBasic() {}
        
        //Initialize a number of filaments
        virtual std::unordered_set<std::unique_ptr<Filament>>* initialize(int numFilaments, int length)
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
            
                Filament* filament = new Filament();
                filament->addSubFilament(FilamentController<NDIM>::_initializer->create(length, maxLength, c_start));
                filaments->emplace(filament);
            }
            return filaments;
        }
        
        ///Extend the front of a filament
        virtual void extendFrontOfFilament(Filament *f)
        {
            ///Find next compartment (1D for now)
            Compartment* cCurrent = f->getFrontSubFilament()->compartment();
            Compartment* cNext = f->getFrontSubFilament()->compartment()->neighbours().back();
            
            ///maxlen, for now
            int maxLength = int(static_cast<CompartmentSpatial<NDIM>*>(cCurrent)->sides()[0] / monomer_size);
            
            SubFilament* s1 = f->getFrontSubFilament();
            SubFilament* s2 = FilamentController<NDIM>::_initializer->create(1, maxLength, cNext);
            
            FilamentController<NDIM>::_initializer->connect(s1,s2);
            ///Add new subfilament
            f->addSubFilament(s2);
        }  
        
    };
    

} //end namespace chem

#endif /* defined(__CytoSim__FilamentController__) */
