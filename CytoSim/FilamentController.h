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
#include "ChemSim.h"

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
        ChemSim &_chem; ///<ptr to chem sim object
        
    public:
        ///Constructor, destructor
        FilamentInitializer(ChemSim &chem) : _chem(chem) {}
        virtual ~FilamentInitializer() {};
        
        ///Set the filament controller
        virtual void setFilamentController (FilamentController<NDIM>* controller)
        { FilamentInitializer<NDIM>::_controller = controller;}
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the SubFilament initialized
        virtual SubFilament* createSubFilament(Filament* parentFilament, int length, int maxlength, Compartment* c) = 0;
        
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
        CompartmentGrid<NDIM>* _grid; ///<compartment grid for updating
        FilamentInitializer<NDIM>* _initializer; ///<initializer, could be any implementation
        
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
    
} //end namespace chem

#endif /* defined(__CytoSim__FilamentController__) */
