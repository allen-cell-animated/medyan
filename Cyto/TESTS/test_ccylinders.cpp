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

TEST(CCylinder, Basic) {
    
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
        
        ccylinder->addInternalReaction(new Reaction<1,1>({sf1, sf2}, 10.0));
    }
    
    ///clone a ccylinder into a new compartment, check its species
    Compartment* newCompartment = new Compartment;
    CCylinder* ccylinderClone = ccylinder->clone(newCompartment);
    
    EXPECT_EQ(20, newCompartment->numberOfSpecies());
    EXPECT_EQ(10, newCompartment->numberOfReactions());
    
    for(int i = 0; i < 10; i++) {
        
        ///check species and reaction equality
        EXPECT_TRUE(ccylinder->getInternalReactions()[i]->is_equal(*ccylinderClone->getInternalReactions()[i]));
        
        EXPECT_TRUE(ccylinderClone->getCMonomer(i)->speciesFilament(0) == ccylinderClone->getCMonomer(i)->speciesFilament(0));
        EXPECT_TRUE(ccylinderClone->getCMonomer(i)->speciesFilament(1) == ccylinderClone->getCMonomer(i)->speciesFilament(1));
    }
    
    ///check destructor
    delete ccylinder;
    EXPECT_EQ(0, c->numberOfSpecies());
    EXPECT_EQ(0, c->numberOfReactions());
}

TEST(CCylinder, AdvancedCloning) {
    
    ChemSim::setInstance(ChemSimInitKey(), new ChemSimpleGillespieImpl);
    ChemSim::initialize(ChemSimInitKey());
    
    Compartment* c1 = new Compartment;
    Compartment* c2 = new Compartment;
    Compartment* c3 = new Compartment;
    
    
    CCylinder* cc1 = new CCylinder(c1);
    CCylinder* cc2 = new CCylinder(c2);
    CCylinder* cc3 = new CCylinder(c3);
    
    ///create cc1 and cc2
    for(int i = 0; i < 10; i++) {
        
        CMonomer* m = new CMonomer;
        SpeciesFilament* sf1 = c1->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("A"));
        SpeciesFilament* sf2 = c1->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("B"));
        
        m->addSpeciesFilament(sf1);
        m->addSpeciesFilament(sf2);
        
        cc1->addCMonomer(m);
        cc1->addInternalReaction(new Reaction<1,1>({sf1, sf2}, 10.0));
    }
    for(int i = 0; i < 10; i++) {
        
        CMonomer* m = new CMonomer;
        SpeciesFilament* sf1 = c2->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("A"));
        SpeciesFilament* sf2 = c2->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName("B"));
        
        m->addSpeciesFilament(sf1);
        m->addSpeciesFilament(sf2);
        
        cc2->addCMonomer(m);
        cc2->addInternalReaction(new Reaction<1,1>({sf1, sf2}, 10.0));
    }
    
    ///add cross-cylinder reactions
    cc1->addCrossCylinderReaction(cc2,
        new Reaction<1,1>({cc1->getCMonomer(9)->speciesFilament(0), cc2->getCMonomer(0)->speciesFilament(0)}, 20.0));
    
    cc2->addCrossCylinderReaction(cc1,
        new Reaction<1,1>({cc2->getCMonomer(0)->speciesFilament(1), cc1->getCMonomer(9)->speciesFilament(1)}, 20.0));
    
    
    EXPECT_EQ(20, c1->numberOfSpecies());
    EXPECT_EQ(11, c1->numberOfReactions());
    
    EXPECT_EQ(20, c2->numberOfSpecies());
    EXPECT_EQ(11, c2->numberOfReactions());
    
    ///now, clone cc1 into c3
    cc3 = cc1->clone(c3);
    EXPECT_EQ(20, c3->numberOfSpecies());
    EXPECT_EQ(11, c3->numberOfReactions());
    
    ///check cross cylinder reaction cloned correctly
    EXPECT_TRUE(cc3->getCrossCylinderReactions()[cc2][0]->is_equal(*cc1->getCrossCylinderReactions()[cc2][0]));
    EXPECT_TRUE(cc2->getCrossCylinderReactions()[cc1][0]->is_equal(*cc2->getCrossCylinderReactions()[cc3][0]));
    
    ///now, delete cc1
    delete cc1;
    
    EXPECT_EQ(0, c1->numberOfSpecies());
    EXPECT_EQ(0, c1->numberOfReactions());
    
    ///reactions involving cc1 should be removed from cc3
    EXPECT_EQ(0, cc3->getCrossCylinderReactions()[cc1].size());
    EXPECT_EQ(cc3->getReactingCylinders().end(),
              std::find(cc3->getReactingCylinders().begin(), cc3->getReactingCylinders().end(), cc1));

}


#endif //DO_THIS_CCYLINDERS_TEST