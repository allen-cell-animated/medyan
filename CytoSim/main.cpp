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
#include "CompartmentContainer.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    
    Compartment *C1 = new Compartment;
    Species *actin = C1->addSpecies("Actin",99U);
    C1->setDiffusionRate(actin,2000);
    Species *profilin = C1->addSpecies("Profilin",29U);
    C1->setDiffusionRate(profilin,2000);
    Species *arp23 = C1->addSpecies("Arp2/3",33U);
    C1->setDiffusionRate(arp23,2000);
    C1->printSpecies();
//    cout << C1->findSpecies(actin_id) << ", " << &C1->findSpecies(actin_id) << endl;
//    cout << C1->species()[actin_id] << ", " << &C1->species()[actin_id] << endl;
    cout << C1->countSpecies() << endl;
    

    vector<Species*> rs1 = {profilin,actin};
    Reaction *r1 = C1->addInternalReaction(rs1, 1, 1, 10.0);
    vector<Species*> rs2 = {actin,profilin};
    Reaction *r2 = C1->addInternalReaction(rs2, 1, 1, 10.0);
    vector<Species*> rs3 = {actin,profilin,arp23};
    Reaction *r3 = C1->addInternalReaction(rs3, 2, 1, 10.0);
    C1->printReactions();
    cout << "Are all Species unique? " << std::boolalpha << C1->areAllSpeciesUnique() << endl;
    
    cout << sizeof(*actin) << ", " << sizeof(actin->getRSpecies()) << ", " << sizeof(*r1) << endl;

    cout << endl << endl;
    
    Compartment *C2 = C1->clone();
    cout << "The clone:" << endl;
    C2->printReactions();
    cout << "End of The Clone" << endl << endl;
    

    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    
    C1->generateDiffusionReactions();
    C2->generateDiffusionReactions();
    
    cout << "C1 Compartment:" << endl;
    C1->printReactions();
    cout << "C1's neigbour size: " << C1->numberOfNeighbours() << endl;
    cout << endl << endl;
    
    cout << "C2 Compartment:" << endl;
    C2->printReactions();
    cout << "C2's neigbour size: " << C2->numberOfNeighbours() << endl;
    cout << endl << endl;
    
    cout << "C3 Compartment:" << endl;
    Compartment C3Val(*C2);
    Compartment *C3 = &C3Val;
    C3->generateDiffusionReactions();
    C3->printReactions();
    cout << "C3's neigbour size: " << C3->numberOfNeighbours() << endl;
    cout << endl << endl;
    
    actin->setN(301);
    cout << "C4 Compartment:" << endl;
    Compartment C4Val;
    C4Val.addSpecies("Ena/VASP",28U);
    C4Val = *C1;
    Compartment *C4 = &C4Val;
    //C3->generateDiffusionReactions();
    C4->printReactions();
    cout << "C4's neigbour size: " << C4->numberOfNeighbours() << endl;
    cout << endl << endl;
    
    SpeciesDiffusing Actin("Actin",12U);
    SpeciesDiffusing Profilin("Profilin",2U);
    Reaction AP({&Actin,&Profilin}, 1, 1, 10.0);
    
    cout << "Finding one Species" << endl;
    Species *C4Actin = C4->findSimilarSpecies(Actin);
    cout << (*C4Actin) << endl;
    
    cout << "Finding one Reaction" << endl;
    Reaction *C4AP = C4->findSimilarInternalReaction(AP);
    cout << (*C4AP) << endl;
    
    cout << "Testing equality operator:" << endl;
    cout << boolalpha << ((*C1)==(*C4)) << endl;
    
    CompartmentCubic<3> CC1;
    vector<float> sides{100.0,100.0,100.0};
    CC1.setSides(sides.begin());
    vector<float> coords{12.3,1.2,22.1};
    CC1.setCoords(coords.begin());
    CC1.printSelf();
    cout << endl << endl;
    
    const int NDIM =3;
    const int NGRID = 2;
    CompartmentsSimpleGrid<NDIM> ccv{NGRID, NGRID, NGRID};
    Compartment &Cproto = ccv.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",9U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",5);
    Cproto.setDiffusionRate(M2,2500);
    vector<Species*> RSM1M2 = {M1,M2};
    Reaction *RM1M2 = Cproto.addInternalReaction(RSM1M2, 1, 1, 311.2);
    ccv.initialize();
    ccv.printSelf();
//    size_t range[NGRID];
//    iota(range,range+NGRID,0);
//    for(size_t i : range)
//        for(size_t j : range)
//            for(size_t k : range)
//            {
//                Compartment *C5 = ccv.getCompartment(i,j,k);
//                cout << i << " " << j << " " << k << ": " << C5 << endl;
//                C5->printSelf();
//                cout << endl;
//            }

    
    cout << "Main exited..." << endl;
    return 0;
}

