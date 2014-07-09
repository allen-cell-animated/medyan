//
//  main.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

/*! \mainpage The Main Page for the CytoSim software package
 
 \section intro_sec Introduction
 
 Cell motility plays a key role in human biology and disease, contributing ubiquitously to such important processes as embryonic development, wound repair and cancer metastasis. Papoian laboratory is interested in gaining deeper understanding of the physical chemistry behind these complex, far-from-equilibrium mechano-chemical processes. His approach is based on combining stochastic reaction-diffusion treatment of cellular biochemical processes with polymer physics of cytoskeletal filament network growth, while explicitly coupling chemistry and mechanics.
 
 Papoian laboratory has developed **CytoSim**, a software package to simulate growth dynamics of actin based filamentous networks *in vitro* and *in vivo*. Recent papers where **CytoSim** or its predecessor, **StochTools**, were used can be found on the publication section of [the Papoian group's main web page: ](http://papoian.chem.umd.edu/)
 
 
 \section install_sec Installation
 
 \subsection step1 Step 1: Prerequisites
 
 The following software packages need to be installed first:
 
 - Boost 1.49
 - GSL ...
 
 \subsection step2 Step 2: Installation of CytoSim itself
 
 Untar the **CytoSim** source code into some directory, enter into the "CytoSim/CytoSim" and execute "make" from the command line.
 
 
 */

#include <iostream>
#include <numeric>
#include <chrono>


//#include "Species.h"
//#include "Reaction.h"
#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"
//#include "ChemSim.h"
//#include "Composite.h"
//#include "SpeciesContainer.h"
#include "FilamentController.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    const int NDIM = 1;
    const int NGRID = 50;
    CompartmentGrid<NDIM> ccv{NGRID};
    CompartmentSpatial<NDIM> &Cproto = ccv.getProtoCompartment();
    
    ///Add diffusing species
    Species *M1 = Cproto.addSpecies("Actin",100U);
    Cproto.setDiffusionRate(M1,2000);

    ///Set side length
    vector<float> sides{100.0};
    Cproto.setSides(sides.begin());
    
    ///init grid
    ccv.initialize();
    
    ///Init filament controller and initializer
    FilamentInitializer<1>* initializer = new ActinBasicInitializer<1>();
    static_cast<ActinBasicInitializer<1>*>(initializer)->setReactionRates();
    
    FilamentController<1>* controller = new FilamentControllerBasic<1>(&ccv, initializer);
    auto filaments = controller->initialize(16, 5);
    
    ///Init chemsim
    ChemNRMImpl chem_sim_impl;
    ChemSim chem(&chem_sim_impl);
    
    ccv.addChemSimReactions(chem);
    chem.initialize();
    
    for(int step = 0; step < 10; step++) {
        
        ///Run 100 steps
        chem.run(100);
        
        ///Print filaments
        for(auto it = filaments->begin(); it != filaments->end(); it++)
            (*it)->printFilament();
    }
    
    
}

