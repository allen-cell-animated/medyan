//
//  CFilamentController.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CFilamentController__
#define __CytoSim__CFilamentController__

#include <iostream>
#include "CFilament.h"
#include "ChemSim.h"

namespace chem {
    
    //Forward declaration
    template<size_t NDIM>
    class CFilamentController;
    
    /// CFilamentInitailizer class is used for initailizing a Csubfilament in a compartment
    /*!
     *  CFilamentInitializer is an abstract class used for initailizing a subfilament, including its
     *  associated species and reactions. SubFilaments can be initialized by constructing
     *  this object and calling its initializer.
     */
    template<size_t NDIM>
    class CFilamentInitializer {
        
    protected:
        CFilamentController<NDIM>* _controller; ///< ptr to controller
        ChemSim &_chem; ///<ptr to chem sim object
        
    public:
        ///Constructor, destructor
        CFilamentInitializer(ChemSim &chem) : _chem(chem) {}
        virtual ~CFilamentInitializer() {};
        
        ///Set the CFilament controller
        virtual void setCFilamentController (CFilamentController<NDIM>* controller)
        { CFilamentInitializer<NDIM>::_controller = controller;}
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the CSubFilament initialized
        ///@param maxlength - length of entire reactive CFilament
        ///@param species - list of species to initialize in CFilament
        virtual CSubFilament* createCSubFilament(CFilament* parentFilament,
                                               Compartment* c,
                                               std::vector<std::string>* species,
                                               int length,
                                               int maxlength) = 0;
        
        ///Connect two CFilaments, back to front
        virtual void connect(CSubFilament* s1, CSubFilament* s2) = 0;
    };
    
    
    /// CFilamentController class is an abstract class for the updating of CFilaments
    /*!
     *  CFilamentController is an abstract class that provides methods to control the updating of
     *  CFilaments after reactions occur. This includes setting up reactions, updating species
     *  and compartments, and creating new Csubfilaments.
     */
    template<size_t NDIM>
    class CFilamentController {

    protected:
        CompartmentGrid<NDIM>* _grid; ///<compartment grid for updating
        CFilamentInitializer<NDIM>* _initializer; ///<initializer, could be any implementation
        
    public:
        
        ///constructor and destructor
        CFilamentController(CompartmentGrid<NDIM>* grid, CFilamentInitializer<NDIM>* initializer)
        : _grid(grid), _initializer(initializer)
        {
            _initializer->setCFilamentController(this);
        }
        
        virtual ~CFilamentController() {}
        
        ///Initialize all CFilaments and reactions, return set of CFilaments
        virtual std::unordered_set<std::unique_ptr<CFilament>>* initialize(int numFilaments, int length) = 0;
        
        ///Extend the front of a CFilament
        virtual void extendFrontOfCFilament(CFilament *f, std::vector<std::string>* species) = 0;
        
//        ///Retract the front of a CFilament
//        virtual void retractFrontOfCFilament(CFilament *f) = 0;
//        
//        ///Retract the back of a CFilament
//        virtual void retractBackOfCFilament(CFilament *f) = 0;
        
    };
    
    ///CFilament REACTION CALLBACKS
    
    template<size_t NDIM>
    struct CFilamentExtensionCallback {
        
        //members
        CFilamentController<NDIM>* _controller;
        CFilament* _filament;
        std::vector<std::string>* _species;
        
        ///Constructor, sets members
        CFilamentExtensionCallback(CFilamentController<NDIM>* controller,
                                  CFilament* filament,
                                  std::vector<std::string>* species) :
            _controller(controller), _filament(filament), _species(species) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            _controller->extendFrontOfCFilament(_filament, _species);
        }
    };
    
} //end namespace chem

#endif /* defined(__CytoSim__CFilamentController__) */
