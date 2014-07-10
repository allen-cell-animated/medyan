//
//  FilamentControllerImpl.h
//  CytoSim
//
//  Created by James Komianos on 7/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__FilamentControllerImpl__
#define __CytoSim__FilamentControllerImpl__

#include <iostream>
#include "FilamentController.h"

namespace chem {
    
    /// FilopodiaInitializer is a basic implementation of the Initializer class that only involves the following:
    
    /// Actin polymerization of (+) end
    /// Actin depolymerization of (+) end
    /// Actin depolymerization of (-) end
    
    /// Capping polymerization of (+) end
    /// Capping depolymerization of (+) end
    
    /// Anti-capping polymerization of (+) end
    /// Anti-capping depolymerization of (+) end
    /// increased rate of actin polymerization with anti-capped (+) end
    
    /// Motor binding to filament
    /// Motor movement (right or left) on filament
    /// Motor loading and unloading
    
    template<size_t NDIM>
    class FilopodiaInitializer : public FilamentInitializer<NDIM> {
        
    private:
        ///REACTION RATES
        float _k_on_plus = 0;
        float _k_off_plus = 0;
        float _k_off_minus = 0;
        
    public:
        
        ///Constructor, does nothing
        FilopodiaInitializer(ChemSim &chem) : FilamentInitializer<NDIM>(chem) {};
        
        ///Destructor, does nothing
        ~FilopodiaInitializer() {};
        
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
        virtual void connect (SubFilament* s1, SubFilament* s2);
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the filament initialized
        ///@param maxlength - length of entire reactive filament
        virtual SubFilament* createSubFilament(Filament* parentFilament,
                                               int length,
                                               int maxlength,
                                               Compartment* c);
        
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
        virtual std::unordered_set<std::unique_ptr<Filament>>* initialize(int numFilaments, int length);
        
        ///Extend the front of a filament
        virtual void extendFrontOfFilament(Filament *f);
        
    };

    
}; //end namespace chem


#endif /* defined(__CytoSim__FilamentControllerImpl__) */
