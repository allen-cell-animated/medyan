
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------


//#define DO_THIS_COMPARTMENT_TEST

#ifdef DO_THIS_COMPARTMENT_TEST

#include "gtest/gtest.h"

#include "common.h"

#include "GController.h"
#include "CompartmentContainer.h"

//Testing neighbor capabilities
TEST(CompartmentTest, Neighbors) {

    Compartment *C1 = new Compartment;
    Compartment *C2 = C1->clone();
    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    
    EXPECT_EQ(1, C1->numberOfNeighbours());
    EXPECT_EQ(1, C2->numberOfNeighbours());
    
    Compartment *C3 = C2->clone();
    C1->addNeighbour(C3);
    
    EXPECT_EQ(2, C1->numberOfNeighbours());
    EXPECT_EQ(0, C3->numberOfNeighbours());
    
    C1->removeNeighbour(C2);
    
    EXPECT_EQ(1, C1->numberOfNeighbours());
}

//Testing species and reaction generation
TEST(CompartmentTest, SpeciesAndReactions){
    
    Compartment *C1 = new Compartment;
    C1->activate();
    
    Species *actin = C1->addSpeciesDiffusing("Actin",99U);
    C1->setDiffusionRate(actin,2000);
    Species *profilin = C1->addSpeciesDiffusing("Profilin",29U);
    C1->setDiffusionRate(profilin,2000);
    Species *arp23 = C1->addSpeciesDiffusing("Arp2/3",33U);
    C1->setDiffusionRate(arp23,2000);
    
    Compartment *C2 = C1->clone();
    
    EXPECT_EQ(2000,C1->getDiffusionRate("Actin"));
    EXPECT_EQ(2000,C2->getDiffusionRate("Actin"));
    
    C1->setDiffusionRate(actin, 9000);
    EXPECT_EQ(9000, C1->getDiffusionRate("Actin"));
    EXPECT_EQ(2000, C2->getDiffusionRate("Actin"));
    
    
    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    C1->generateAllDiffusionReactions();
    C2->generateAllDiffusionReactions();
    
    EXPECT_EQ(3, C1->numberOfSpecies());
    EXPECT_EQ(3, C2->numberOfSpecies());
    
    EXPECT_EQ(0, C1->numberOfInternalReactions());
    EXPECT_EQ(3, C1->numberOfReactions());
    EXPECT_EQ(0, C2->numberOfInternalReactions());
    EXPECT_EQ(3, C2->numberOfReactions());

}

//Testing neighbors, species and reaciton generation
TEST(CompartmentContainerTest, Main) {
    
    
    SystemParameters::GParams.compartmentSizeX = 100.0;
    SystemParameters::GParams.compartmentSizeY = 100.0;
    SystemParameters::GParams.compartmentSizeZ = 100.0;
    
    SystemParameters::GParams.NX = 50;
    SystemParameters::GParams.NY = 50;
    SystemParameters::GParams.NZ = 50;
    
    SystemParameters::GParams.monomerSize = 2.7;
    SystemParameters::GParams.cylinderSize = 27.0;
    
    SystemParameters::GParams.nDim = 3;
    
    int _numSpecies; ///for testing

    GController g;
    g.initializeGrid();
    
    Compartment &Cproto = CompartmentGrid::instance()->getProtoCompartment();
    
    Species *M1 = Cproto.addSpeciesDiffusing("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpeciesDiffusing("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    
    Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
    
    ///initialize all compartments with species
    for(auto &c : CompartmentGrid::instance()->children())
    {
        Compartment *C = static_cast<Compartment*>(c.get());
        *C = Cproto;
    }
    
    ///activate all compartments for diffusion
    CompartmentGrid::instance()->activateAll();
    
    ///Generate all diffusion reactions
    for(auto &c : CompartmentGrid::instance()->children())
    {
        Compartment *C = static_cast<Compartment*>(c.get());
        C->generateAllDiffusionReactions();
    }
    
    _numSpecies = CompartmentGrid::instance()->countSpecies();
    
    Compartment *C1 = GController::getCompartment(vector<size_t>{10U,10U,10U});
    
    EXPECT_EQ(6, C1->numberOfNeighbours());
    EXPECT_EQ(2, C1->numberOfSpecies());
    EXPECT_EQ(14, C1->numberOfReactions());
    
    Compartment*C2 = GController::getCompartment(vector<size_t>{10U,11U,10U});
    C1->removeNeighbour(C2);
    
    EXPECT_EQ(5, C1->numberOfNeighbours());
    
}

#endif // DO_THIS_COMPARTMENT_TEST






