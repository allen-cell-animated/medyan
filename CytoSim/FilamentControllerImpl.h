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
    
    /// formin polymerization of (+) end
    /// formin depolymerization of (+) end
    /// increased rate of actin polymerization with formin (+) end
    
    /// Motor binding to filament
    /// Motor movement (right or left) on filament
    /// Motor loading and unloading
    
    template<size_t NDIM>
    class FilopodiaInitializer : public FilamentInitializer<NDIM> {
        
    private:
        ///REACTION RATES
        //basic
        float _k_on_plus = 21.0;
        float _k_off_plus = 1.4;
        float _k_off_minus = 1.4;
        
        //capping
        float _k_capping_on_plus = 50.0;
        float _k_capping_off_plus = 0.06;
        
        //formin
        float _k_formin_on_plus = 10.0;
        float _k_formin_off_plus = 1.4;
        float _k_accel_on_plus = 100.0;
        
        //motors
        float _k_binding = 19.0;
        float _k_unbinding = 10.0;
        float _k_forward_step = 50.0;
        float _k_backward_step = 5.0;
        float _k_load = 20.0;
        float _k_unload = 10.0;
        
        
    public:
        
        ///Constructor, does nothing
        FilopodiaInitializer(ChemSim &chem) : FilamentInitializer<NDIM>(chem) {};
        
        ///Destructor, does nothing
        ~FilopodiaInitializer() {};
        
//        ///Set reaction rates for this subfilament
//        virtual void setReactionRates(float kOnPlus = 21.0,
//                                      float kOffPlus = 1.4,
//                                      float kOffMinus = 1.4
//                                      float kCappingOnPlus = 50.0
//                                      float kCappingOffPlus = 0.06
//                                      float kForminOnPlus = 10.0
//                                      float kForminOffPlus = 1.4
//                                      float kAccelOnPlus = 100.0)
//        {
//            _k_on_plus = kOnPlus;
//            _k_off_plus = kOffPlus;
//            _k_off_minus = kOffMinus;
//            _k_capping_on_plus = kCappingOnPlus;
//            _k_capping_off_plus = kCappingOffPlus;
//            _k_formin_on_plus = kForminOnPlus;
//            _k_formin_off_plus = 0;
//            _k_accel_on_plus = 0;
//            
//        }
        
        ///Connect two filaments, back to front
        ///For this impl, only add a polymerization reaction between them
        virtual void connect (SubFilament* s1, SubFilament* s2);
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the filament initialized
        ///@param maxlength - length of entire reactive filament
        ///@param species - list of species to initialize in filament
        virtual SubFilament* createSubFilament(Filament* parentFilament,
                                               Compartment* c,
                                               std::vector<std::string>* species,
                                               int length,
                                               int maxlength);
        
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
        virtual void extendFrontOfFilament(Filament *f, std::vector<std::string>* species);
        
    };

    
}; //end namespace chem


#endif /* defined(__CytoSim__FilamentControllerImpl__) */
