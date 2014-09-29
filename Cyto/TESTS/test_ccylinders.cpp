//
//  test_ccylinders.cpp
//  Cyto
//
//  Created by James Komianos on 9/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
#define DO_THIS_CCYLINDERS_TEST

#ifdef DO_THIS_CCYLINDERS_TEST

#include <iostream>
#include "gtest/gtest.h"

#include "common.h"
#include "CCylinder.h"

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
    
    
    
    
    
}








#endif //DO_THIS_CCYLINDERS_TEST