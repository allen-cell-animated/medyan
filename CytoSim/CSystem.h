//
//  CSystem.h
//  CytoSim
//
//  Created by James Komianos on 7/21/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CSystem__
#define __CytoSim__CSystem__

#include <iostream>
#include "CompartmentContainer.h"
#include "CSystem.h"


namespace chem {

    template <size_t NDIM>
    class CSystem {
        
//    private:
//        CompartmentGrid<NDIM>* _grid;
//        CSystem<NDIM>* _controller;
//        
//    public:
//        CSystem(CSystem<NDIM> *controller) : _controller(controller)
//        {
//            
//            _controller->setSystem(this);
//        }
//        
//        ~CSystem()
//        {
//            delete _controller;
//            delete _grid;
//        }
//        
//        void getCSystem() {return _controller;}
//        void getCompartmentGrid() {return _grid;}
//        
//        
        
    protected:
        CFilamentInitializer<NDIM>* _initializer; ///<initializer, could be any implementation
        CompartmentGrid<NDIM>* _grid;
        
        std::unordered_set<std::unique_ptr<CFilament>> _filaments;///< filaments that this is controlling
        
    public:
        
        ///constructor and destructor
        CSystem() {}
        
        ///delete filaments
        virtual ~CSystem(){}
        
        virtual setInitializer(CFilamentInitializer<NDIM>* initializer)
        {
            _initializer = initializer;
        }
        
        virtual initializeGrid(std::initializer_list<size_t> gridSize)
        {
            _grid = new CompartmentGrid<NDIM>(gridSize);
        }
        
        
        
        
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
    
};


#endif /* defined(__CytoSim__CSystem__) */
