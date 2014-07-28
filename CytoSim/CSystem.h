//
//  CSystem.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CSystem__
#define __CytoSim__CSystem__

#include <iostream>
#include "CFilament.h"
#include "ChemSim.h"

namespace chem {
    
    //Forward declaration
    template<size_t NDIM>
    class CSystem;
    
    /// CFilamentInitailizer class is used for initailizing a Csubfilament in a compartment
    /*!
     *  CFilamentInitializer is an abstract class used for initailizing a subfilament, including its
     *  associated species and reactions. SubFilaments can be initialized by constructing
     *  this object and calling its initializer.
     */
    template<size_t NDIM>
    class CFilamentInitializer {
        
    protected:
        CSystem<NDIM>* _csystem; ///< ptr to controller
        ChemSim &_chem; ///<ptr to chem sim object
        
    public:
        ///Constructor, destructor
        CFilamentInitializer(ChemSim &chem) : _chem(chem) {}
        virtual ~CFilamentInitializer() {};
        
        ///Set the CFilament controller
        virtual void setCSystem (CSystem<NDIM>* csystem) {_csystem = csystem;}
        
        ///Initialize proto compartment based on this implementation
        virtual void initializeProtoCompartment(CompartmentSpatial<NDIM>& Cproto) = 0;
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the CSubFilament initialized
        ///@param maxlength - length of entire reactive CFilament
        ///@param species - list of species to initialize in CFilament
        virtual CSubFilament* createCSubFilament(CFilament* parentFilament,
                                               Compartment* c,
                                               std::vector<std::string> species,
                                               int length) = 0;
        
        ///Remove a CSubFilament, based on the given simulation
        virtual void removeCSubFilament(CFilament* parentFilament) = 0;
        
        ///Update based on a given reaction occuring
        virtual void update(CFilament* f, ReactionBase* r) = 0;
        
        ///get this chemsim engine
        ChemSim& getChemSim() {return _chem;}
    };
    
    
    /// CSystem class is an abstract class for the updating of CFilaments
    /*!
     *  CSystem is an abstract class that provides methods to control the updating of
     *  CFilaments after reactions occur. This includes setting up reactions, updating species
     *  and compartments, and creating new Csubfilaments.
     */
    template<size_t NDIM>
    class CSystem {

    protected:
        CompartmentGrid<NDIM>* _grid; ///<compartment grid
        CFilamentInitializer<NDIM>* _initializer; ///<initializer, could be any implementation
        
        std::vector<std::unique_ptr<CFilament>> _filaments;///< filaments that this is controlling
        
    public:
        
        ///constructor and destructor
        CSystem(CFilamentInitializer<NDIM>* initializer) : _initializer(initializer)
        {
            
            ///Set up grid
            if(NDIM ==1)
                _grid = new CompartmentGrid<NDIM>({NGRID});
            
            else if(NDIM == 2)
                _grid = new CompartmentGrid<NDIM>({NGRID,NGRID});
            
            else
                _grid = new CompartmentGrid<NDIM>({NGRID,NGRID,NGRID});
            
            ///Set system
            _initializer->setCSystem(this);
            
            ///init protocompartment
            CompartmentSpatial<NDIM> &Cproto = _grid->getProtoCompartment();
            _initializer->initializeProtoCompartment(Cproto);
            
            ///init grid
            _grid->initialize();
        }
        
        ///Init chemistry in grid (including adding diffusion reactions, updating filament reactions)
        void initChem()
        {
            ///Generate diffusion reactions
            _grid->generateDiffusionReactions();
            
            ///Init chemsim
            ChemSim& chem = _initializer->getChemSim();
            _grid->addChemSimReactions(chem);
            
            chem.initialize();
            
            ///loop through filaments, passivate/activate reactions
            for (auto &f : _filaments)
                f->updateReactions();
            
            //chem.printReactions();
        }
        
        ///delete filaments
        virtual ~CSystem(){}
        
        ///Initialize a CFilament
        virtual CFilament* initializeCFilament(int length) = 0;
        
        ///Extend the front of a CFilament
        virtual void extendFrontOfCFilament(CFilament *f, std::vector<std::string> species) = 0;
        
        ///Retract the front of a CFilament
        virtual void retractFrontOfCFilament(CFilament *f) = 0;
        
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
        
        ///get grid
        virtual CompartmentGrid<NDIM>* grid() {return _grid;}
        
    };
    
    ///REACTION CALLBACKS
    
    ///Extension callback
    template<size_t NDIM>
    struct CFilamentExtensionCallback {
        
        //members
        CSystem<NDIM>* _csystem;
        CFilament* _filament;
        std::vector<std::string> _species;
        
        ///Constructor, sets members
        CFilamentExtensionCallback(CSystem<NDIM>* csystem,
                                   CFilament* filament,
                                   std::vector<std::string> species) :
        _csystem(csystem), _filament(filament), _species(species) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            _csystem->extendFrontOfCFilament(_filament, _species);
            _filament->getFrontCSubFilament()->updateReactions();
            _csystem->update(_filament, r);
        }
    };
    
    ///Retraction callback
    template<size_t NDIM>
    struct CFilamentRetractionCallback {
        
        //members
        CSystem<NDIM>* _csystem;
        CFilament* _filament;
        
        ///Constructor, sets members
        CFilamentRetractionCallback(CSystem<NDIM>* csystem,
                                   CFilament* filament) :
        _csystem(csystem), _filament(filament) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            _csystem->retractFrontOfCFilament(_filament);
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
    
    
    
} //end namespace chem

#endif /* defined(__CytoSim__CSystem__) */
