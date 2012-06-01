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
#include <fstream>

#include <boost/heap/pairing_heap.hpp>

#include "System.h"
#include "ChemNRMImpl.h"
#include "ChemSim.h"
#include "Signaling.h"

#include "Composite.h"

using namespace std;
using namespace chem;


int main(int argc, const char * argv[])
{
    auto C = make_unique<Composite>();
    cout << C->getFullName() << endl;
    C->addSpeciesUniq(make_unique<SpeciesBulk>("G-Actin",42));
    C->addSpecies<SpeciesBulk>("Profilin",31);
    
    auto a1 = C->species(0);
    auto a2 = C->species(1);
    
    cout << *a1 << "\n" << *a2 << endl;
    
    for(auto &s: C->species()){
        cout << s->getFullName() << s->getN() << endl; 
        cout << &s << endl;
    }
    
    cout << "\n\n\n";
    
    auto D = make_unique<Composite>();
    C->addChild(std::move(D));
    
    Composite *F = C->getRoot();
    F->addSpecies<SpeciesBulk>("Arp2/3",22);

    
    for(auto &s: F->species()){
        cout << s->getFullName() << s->getN() << endl; 
        cout << &s << endl;
    }
    
    cout << "F has " << F->countSpecies() << " species" << endl;
    
    cout << "Main exited..." << endl;
    return 0;
}

