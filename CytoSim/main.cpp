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

//#include "Species.h"
//#include "Reaction.h"
//#include "ChemNRMImpl.h"
////#include "ChemGillespieImpl.h"
////#include "ChemSimpleGillespieImpl.h"
//#include "ChemSim.h"

//#include "Composite.h"

//#include "SpeciesContainer.h"

#include "Compartment.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    
    Compartment C;
    size_t actin_id = C.addSpecies("Actin",species_copy_t(99));
    size_t profilin_id = C.addSpecies("Profilin",species_copy_t(29));
    size_t arp23_id = C.addSpecies("Arp2/3",species_copy_t(33));
    C.printSpecies();
//    cout << C.findSpecies(actin_id) << ", " << &C.findSpecies(actin_id) << endl;
//    cout << C.species()[actin_id] << ", " << &C.species()[actin_id] << endl;
    cout << C.countSpecies() << endl;
    
    Species &S(C.species()[actin_id]);
    
    cout << sizeof(S);
    
    cout << "Main exited..." << endl;
    return 0;
}

