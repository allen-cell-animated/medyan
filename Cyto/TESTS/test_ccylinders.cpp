//
//  test_ccylinders.cpp
//  Cyto
//
//  Created by James Komianos on 9/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
//#define DO_THIS_CCYLINDERS_TEST

#ifdef DO_THIS_CCYLINDERS_TEST

#include <iostream>
#include "gtest/gtest.h"

#include "common.h"
#include "CCylinder.h"
#include "ChemSimpleGillespieImpl.h"

TEST(CMonomer, Main) {
    
    CMonomer* m1 = new CMonomer;
    Compartment* c1 = new Compartment;
    Compartment* c2 = new Compartment;
    
    SpeciesFilament* sf1 = c1->addSpeciesFilament("A");
    SpeciesFilament* sf2 = c1->addSpeciesFilament("B");
    
    m1->addSpeciesFilament(sf1);
    m1->addSpeciesFilament(sf2);
    
    ///copy constructor
    CMonomer* m2 = new CMonomer(*m1, c2);
    
    Species* sf3 = c2->findSimilarSpecies(*sf1);
    Species* sf4 = c2->findSimilarSpecies(*sf2);
    
    EXPECT_EQ(sf3, m2->speciesFilamentVector()[0]);
    EXPECT_EQ(sf4, m2->speciesFilamentVector()[1]);
    
    EXPECT_EQ(sf3, m2->speciesFilament(0));
    EXPECT_EQ(sf4, m2->speciesFilament(1));
}

TEST(CCylinder, Main) {
    
    ChemSim::setInstance(ChemSimInitKey(), new ChemSimpleGillespieImpl);
    ChemSim::initialize(ChemSimInitKey());
    
    Compartment* c = new Compartment;
    SpeciesFilament* sf1;
    SpeciesFilament* sf2;
    CCylinder* ccylinder = new CCylinder(c);
    
    for(int i = 0; i < 10; i++) {
    
        CMonomer* m = new CMonomer;
        sf1 = c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("A"));
        sf2 = c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("B"));
        
        m->addSpeciesFilament(sf1);
        m->addSpeciesFilament(sf2);
        
        ccylinder->addCMonomer(m);
        
        ccylinder->addReaction(new Reaction<1,1>({sf1, sf2}, 10.0));
    }
    
    ///clone a ccylinder into a new compartment, check its species
    Compartment* newCompartment = new Compartment;
    CCylinder* ccylinderClone = ccylinder->clone(newCompartment);
    
    EXPECT_EQ(20, newCompartment->numberOfSpecies());
    EXPECT_EQ(10, newCompartment->numberOfReactions());
    
    for(int i = 0; i < 10; i++) {
        EXPECT_TRUE(ccylinderClone->getCMonomer(i)->speciesFilament(0) == ccylinderClone->getCMonomer(i)->speciesFilament(0));
        EXPECT_TRUE(ccylinderClone->getCMonomer(i)->speciesFilament(1) == ccylinderClone->getCMonomer(i)->speciesFilament(1));
    }
    
    ///check destructor
    delete ccylinder;
    EXPECT_EQ(0, c->numberOfSpecies());
    EXPECT_EQ(0, c->numberOfReactions());
}


#endif //DO_THIS_CCYLINDERS_TEST