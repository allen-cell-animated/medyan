//
//  ChemInitializerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemInitializerImpl__
#define __Cyto__ChemInitializerImpl__

#include <iostream>
#include <vector>

namespace chem {
    class Compartment;
    class CCylinder;
    
    ///ChemInitializerImpl is an abstract base class for initialization of all chemistry in the system
    
    /*  
     *  Specific initializers should inherit from ChemInitializerImpl. A user will then attach the corresponding 
     *  initializer to ChemInitializer via the initializer base class, ChemInitializerImpl.
     */
    public ChemInitializerImpl {
        
    public:
        ///Initializer, based on the given simulation
        ///@param length - starting length of the CCylinder initialized
        ///@param species - list of species to initialize in CCylinder
        virtual CCylinder* createCCylinder(Compartment* c, std::vector<std::string> species, int length) = 0;
        
        ///Remove a CCylinder, based on the given simulation
        virtual void removeCCylinder(CCylinder *cylinder) = 0;

    };
    
}; //end namespace chem

#endif /* defined(__Cyto__ChemInitializerImpl__) */
