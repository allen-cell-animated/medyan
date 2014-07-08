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
    
    
    ReactionBase *r1 = C1->addInternal<Reaction,1,1>({profilin,actin}, 10.0);
    C1->addInternal<Reaction,1,1>({actin,profilin}, 10.0);
    C1->addInternalReaction<2,1>(vector<Species*>{actin,profilin,arp23}, 10.0);
    
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
    Reaction<1,1> AP({&Actin,&Profilin}, 10.0);
    
    cout << "Finding one Species" << endl;
    Species *C4Actin = C4->findSimilarSpecies(Actin);
    cout << (*C4Actin) << endl;
    
    cout << "Finding one Reaction" << endl;
    ReactionBase *C4AP = C4->findSimilarInternalReaction(AP);
    cout << (*C4AP) << endl;
    
    cout << "Testing equality operator:" << endl;
    cout << boolalpha << ((*C1)==(*C4)) << endl;
    
    CompartmentSpatial<3> CC1;
    vector<float> sides{100.0,100.0,100.0};
    CC1.setSides(sides.begin());
    vector<float> coords{12.3,1.2,22.1};
    CC1.setCoords(coords.begin());
    CC1.printSelf();
    cout << endl << endl;
    
    auto chk_ccv_1 = chrono::high_resolution_clock::now();
    const int NDIM =3;
    const int NGRID = 50;
    CompartmentGrid<NDIM> ccv{NGRID, NGRID, NGRID};
    CompartmentSpatial<NDIM> &Cproto = ccv.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    Cproto.setSides(sides.begin());
    
    ReactionBase *RM1M2 = Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
    Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    ccv.initialize();
    
    //ccv.printSelf();
    cout << "Num of Species: " << ccv.countSpecies() << endl;
    cout << "Num of Reactions: " << ccv.countReactions() << endl;
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
    
    
    //ChemSimpleGillespieImpl chem_sim_impl;
    //ChemGillespieImpl chem_sim_impl;
    ChemNRMImpl chem_sim_impl;
    ChemSim chem(&chem_sim_impl);
    
    cout << "ccv.addChemSimReactions(chem)..." << endl;
    ccv.addChemSimReactions(chem);
    cout << "chem.initialize()..." << endl;
    chem.initialize();
    Compartment *C000 = ccv.getCompartment(0U,0U,0U);
    C000->printSelf();
    
    //chem.printReactions();
    long int NUM_OF_STEP = pow(10,7);
    //chrono::high_resolution_clock::time_point chk1, chk2;
    auto chk1 = chrono::high_resolution_clock::now();
    chem.run(NUM_OF_STEP);
    auto chk2 = chrono::high_resolution_clock::now();
    
    
    C000->printSelf();
    cout << "tau=" << tau() <<endl;
    
    chrono::duration<double> elapsed_run(chk2-chk1); // duration in seconds
    //long long int = std::chrono::duration_cast<std::chrono::nanoseconds>(chk2-chk1).count();
    
    cout << "Main Elapsed for ChemNRMImpl::run(...): dt=" << elapsed_run.count() << endl;
    
    auto chk_ccv_2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_ccv(chk_ccv_2-chk_ccv_1); // duration in seconds
    cout << "Main Total Elapsed for CompartmentGrid<NDIM>...: dt=" << elapsed_ccv.count() << endl;
    
    cout << "Size of RSpecies: " << sizeof(M1->getRSpecies()) << endl;
    cout << "Size of ReactionBase: " << sizeof(*RM1M2) << endl;
    cout << "Size of Reaction: " << sizeof(*static_cast<Reaction<1,1>*>(RM1M2)) << endl;
    
    cout << "Main exited..." << endl;
    
    return 0;
}

