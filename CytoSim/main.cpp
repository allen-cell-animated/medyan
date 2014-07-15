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
#include "CFilopodiaController.h"
#include "CMembrane.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    const int NDIM = 1;
    const int NGRID = 5;
    CompartmentGrid<NDIM> ccv{NGRID};
    CompartmentSpatial<NDIM> &Cproto = ccv.getProtoCompartment();
    
    ///Add diffusing species
    Species *M1 = Cproto.addSpecies("Actin",10U);
    Species *M2 = Cproto.addSpecies("Capping",0U);
    Species *M3 = Cproto.addSpecies("X-Formin",5U);
    Species *M4 = Cproto.addSpecies("Myosin", 5U);
    Cproto.setDiffusionRate(M1,2000);
    Cproto.setDiffusionRate(M2,2000);
    Cproto.setDiffusionRate(M3,2000);
    Cproto.setDiffusionRate(M4,2000);

    ///Set side length
    vector<float> sides{100.0};
    Cproto.setSides(sides.begin());
    
    ///init grid
    ccv.initialize();
    
    ///Create chemsim and init
    ChemNRMImpl chem_sim_impl;
    ChemSim chem(&chem_sim_impl);
    ccv.addChemSimReactions(chem);
    chem.initialize();
    
    CMembrane mem;
    
    ///Init filament initializer
    FilopodiaInitializer<1> initializer{chem, mem};
    
    //Init filament controller
    CFilamentControllerFilopodia<1> controller{&ccv, &initializer};
    controller.initialize(16, 20);
    
    //chem.printReactions();
    
    for(int step = 0; step < 1; step++) {
        
        ///Run 100 steps
        chem.run(10000);
        
        controller.printFilaments();
       
    }
    
    std::cout << "Done!" <<std::endl;
    
}

