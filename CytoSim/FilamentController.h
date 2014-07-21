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
#include "CFilament.h"
#include "ChemSim.h"

namespace chem {
    
    //Forward declaration
    template<size_t NDIM>
    class FilamentController;
    
    /// CFilamentInitailizer class is used for initailizing a Csubfilament in a compartment
    /*!
     *  CFilamentInitializer is an abstract class used for initailizing a subfilament, including its
     *  associated species and reactions. SubFilaments can be initialized by constructing
     *  this object and calling its initializer.
     */
    template<size_t NDIM>
    class CFilamentInitializer {
        
    protected:
        FilamentController<NDIM>* _controller; ///< ptr to controller
        ChemSim &_chem; ///<ptr to chem sim object
        
    public:
        ///Constructor, destructor
        CFilamentInitializer(ChemSim &chem) : _chem(chem) {}
        virtual ~CFilamentInitializer() {};
        
        ///Set the CFilament controller
        virtual void setFilamentController (FilamentController<NDIM>* controller)
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
        
        //Update based on a given reaction occuring
        virtual void update(CFilament* f, ReactionBase* r) = 0;
    };
    
    
    /// FilamentController class is an abstract class for the updating of CFilaments
    /*!
     *  FilamentController is an abstract class that provides methods to control the updating of
     *  CFilaments after reactions occur. This includes setting up reactions, updating species
     *  and compartments, and creating new Csubfilaments.
     */
    template<size_t NDIM>
    class FilamentController {

    protected:
        CompartmentGrid<NDIM>* _grid; ///<compartment grid for updating
        CFilamentInitializer<NDIM>* _initializer; ///<initializer, could be any implementation
        
        std::unordered_set<std::unique_ptr<CFilament>> _filaments;///< filaments that this is controlling
        
    public:
        
        ///constructor and destructor
        FilamentController(CompartmentGrid<NDIM>* grid, CFilamentInitializer<NDIM>* initializer)
        : _grid(grid), _initializer(initializer)
        {
            _initializer->setFilamentController(this);
        }
        
        ///delete filaments
        virtual ~FilamentController(){}
        
        ///Initialize all CFilaments and reactions, return set of CFilaments
        virtual void initialize(int numFilaments, int length) = 0;
        
        ///Extend the front of a CFilament
        virtual void extendFrontOfCFilament(CFilament *f, std::vector<std::string>* species) = 0;
        
        ///Update based on a given reaction occuring
        void update(CFilament* f, ReactionBase* r) {
            _initializer->update(f, r);
        }
        
        ///Print filaments
        virtual void printFilaments() {
            int index = 0;
            for(auto it = _filaments.begin(); it != _filaments.end(); it++) {
                std::cout << "FILAMENT " << index++ << std::endl;
                (*it)->printCFilament();
                std::cout << std::endl;
            }
        }
        
        
        
    };
    
} //end namespace chem

#endif /* defined(__CytoSim__FilamentController__) */
