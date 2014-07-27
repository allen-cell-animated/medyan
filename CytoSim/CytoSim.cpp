//
//  CytoSim.cpp
//  CytoSim
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

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"

#include "FilopodiaCSystem.h"
#include "CMembrane.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    ///Create chemsim and init
    ChemNRMImpl chem_sim_impl;
    ChemSim chem(&chem_sim_impl);

    CMembrane mem;
    
    ///Init filament initializer
    SimpleInitializer<1> initializer{chem, mem};
    
    //Init filament controller
    FilopodiaCSystem<1> csystem{&initializer};
    
    
    for(int i = 0; i < 16; i++)
        csystem.initializeCFilament(81.0);
    
    //csystem.printFilaments();
    
    //chem.printReactions();
    
    for(int step = 0; step < 100; step++) {
        
        ///Run 100 steps
        auto chk1 = chrono::high_resolution_clock::now();
        chem.run(80000);
        auto chk2 = chrono::high_resolution_clock::now();
        
        csystem.retrogradeFlow();
        
        chrono::duration<double> elapsed_run(chk2-chk1); // duration in seconds
        //long long int = std::chrono::duration_cast<std::chrono::nanoseconds>(chk2-chk1).count();
       
    }
    csystem.printFilaments();
    cout << "tau=" << tau() <<endl;
    std::cout << "Done!" <<std::endl;
    
}

