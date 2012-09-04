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
#include "ReactionContainer.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    
    Compartment C;
    Species *actin = C.addSpecies("Actin",99U);
    Species *profilin = C.addSpecies("Profilin",29U);
    Species *arp23 = C.addSpecies("Arp2/3",33U);
    C.printSpecies();
//    cout << C.findSpecies(actin_id) << ", " << &C.findSpecies(actin_id) << endl;
//    cout << C.species()[actin_id] << ", " << &C.species()[actin_id] << endl;
    cout << C.countSpecies() << endl;
    
    vector<Species*> rs1 = {actin,profilin};
    Reaction *r1 = C.addReaction(rs1, 1, 1, 10.0);
    vector<Species*> rs2 = {profilin,actin};
    Reaction *r2 = C.addReaction(rs2, 1, 1, 10.0);
    vector<Species*> rs3 = {actin,profilin,arp23};
    Reaction *r3 = C.addReaction(rs3, 2, 1, 10.0);
    C.printReactions();
    cout << "Are all Species unique? " << std::boolalpha << C.areAllSpeciesUnique() << endl;
    
    cout << sizeof(*actin) << ", " << sizeof(actin->getRSpecies()) << ", " << sizeof(*r1) << endl;

    cout << endl << endl;
    
    Compartment *C2 = C.clone();
    cout << "The clone:" << endl;
    C2->printReactions();

    
    
    cout << "Main exited..." << endl;
    return 0;
}

