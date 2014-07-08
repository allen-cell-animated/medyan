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
#include "CompartmentContainer.h"
#include "Filament.h"

namespace chem {
    
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
        
    public:
        ///Initialize all filaments and reactions, return set of filaments
        virtual std::unordered_set<std::unique_ptr<Filament>>& initialize(int numFilaments, int length) = 0;
        
        ///Extend the front of a filament
        virtual void extendFrontOfFilament(Filament *f) = 0;
        
        ///Retract the front of a filament
        virtual void retractFrontOfFilament(Filament *f) = 0;
        
        ///Retract the back of a filament
        virtual void retractBackOfFilament(Filament *f) = 0;
        
    };
    
    
    ///FILAMENT REACTION CALLBACKS
    
    
    struct FilamentReactionCallback {
        
        FilamentReactionCallback(FilamentController<1>* controller, Filament* f)
        : _controller(controller), _filament_to_update(f) {}
        
        void operator() (ReactionBase *r){
            _controller->extendFrontOfFilament(_filament_to_update);
            //        cout << "ReactionCallback was called by\n" << *r << endl;
        }
        Filament* _filament_to_update;
        FilamentController<1>* _controller;
    };
    

    /// FilamentControllerBasic is a basic implementation for updating filaments
    template <size_t NDIM>
    class FilamentControllerBasic : public FilamentController<NDIM> {
      
        ///Default constructor
        FilamentControllerBasic(CompartmentGrid<NDIM>* grid) : FilamentController<NDIM>::_grid(grid) {}
        
        ///Default destructor, does nothing
        ~FilamentControllerBasic() {}
        
        //Initialize
        virtual std::unordered_set<std::unique_ptr<Filament>>& initalize(int numFilaments, int length)
        {
            Compartment* c_start;
            ///Starting compartment for 1D, all filaments start in compartment 0
            if (NDIM == 1) {
                c_start = FilamentController<NDIM>::_grid->getCompartment(0);
            }
            else {
                std::cout << "Multiple dimensional implementation not done. Exiting." << std::endl;
                exit(EXIT_FAILURE);
            }
            
            //Look up diffusing species in this compartment
            Species* actinDiffusing = c_start->findSpeciesByName("Actin");
            Species* forminDiffusing = c_start->findSpeciesByName("Formin");
            Species* cappingDiffusing = c_start->findSpeciesByName("Capping");
            
            
            std::unordered_set<std::unique_ptr<Filament>> _filaments;
            for (int i = 0; i < numFilaments; i++)
            {
                ///Create filament
                Filament *f = new Filament(c_start);
                SubFilament* subf = f->getFrontSubFilament();
                
                ///Add monomer
                for (int j = 0; j < length; j++)
                {
                    ///Add species we want
                    Species *actinFilament = c_start->addSpeciesFilament("Actin", 1);
                    Species *forminFilament = c_start->addSpeciesFilament("Formin", 0);
                    Species *cappingFilament = c_start->addSpeciesFilament("Capping", 0);
                    
                    ///Front and back
                    Species *back = c_start->addSpeciesFilament("Back", 0);
                    Species *front = c_start->addSpeciesFilament("Front", 0);
                    if(j == 0) back->setN(1);
                    if(j == length - 1) front->setN(1);
                        
                    ///add to subfilament
                    subf->addMonomer(new Monomer({actinFilament, forminFilament,
                                            cappingFilament, back, front}, c_start));
                }
                
                ///Add polymerization reactions
                for (int j = 0; j < length; j++) {
                    
                    if (j != length - 1) {
                        c_start->addInternal<Reaction,2,2>({subf->monomer(j).species("Front"), actinDiffusing,
                            subf->monomer(j+1).species("Actin"), subf->monomer(j+1).species("Front")}, 10.0);
                    }
                    else {
                        
                        FilamentReactionCallback* callback = new FilamentReactionCallback(this, f);
                        
                        ReactionBase* r =
                        c_start->addInternal<Reaction,2,0>({subf->monomer(j).species("Front"), actinDiffusing}, 10.0);
                        boost::signals2::shared_connection_block rcb(r->connect((*callback)()), false);
                    }
                        
                }
                
                
            }
            
            
            
        }
        
        ///Extend the front of a filament
        virtual void extendFrontOfFilament(Filament *f) = 0;
        
        
        
        
    };
    

} //end namespace chem

#endif /* defined(__CytoSim__FilamentController__) */
