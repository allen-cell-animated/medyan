//
//  test_boundaries.cpp
//  Cyto
//
//  Created by James Komianos on 9/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

//#define DO_THIS_BOUNDARY_TEST

#ifdef DO_THIS_BOUNDARY_TEST

#include <iostream>
#include "gtest/gtest.h"

#include "BoundaryImpl.h"
#include "SystemParameters.h"

TEST(BoundaryElementTest, Main) {
    
    SystemParameters::GParams.compartmentSizeX = 10.0;
    SystemParameters::GParams.compartmentSizeY = 10.0;
    SystemParameters::GParams.compartmentSizeZ = 10.0;
    
    SystemParameters::GParams.NX = 5;
    SystemParameters::GParams.NY = 5;
    SystemParameters::GParams.NZ = 5;
    
    SystemParameters::GParams.nDim = 3;
    
    GController::initializeGrid();
    
    BoundaryElement* b1 = new PlaneBoundaryElement({15.0,15.0,15.0}, {1,0,0}, 1.0);
    
    ///compartment
    EXPECT_EQ(GController::getCompartment(std::vector<double>{15.0,15.0,15.0}), b1->getCompartment());
    
    ///adding and removing neighbors
    BoundaryElement* b2 = new PlaneBoundaryElement({25.0,15.0,15.0}, {1,0,0}, 1.0);
    BoundaryElement* b3 = new PlaneBoundaryElement({35.0,15.0,15.0}, {1,0,0}, 1.0);
    
    b1->addNeighbor(b2);
    b2->addNeighbor(b1);
    b1->addNeighbor(b3);
    
    EXPECT_TRUE(b1->isNeighbor(b2));
    EXPECT_TRUE(b2->isNeighbor(b1));
    EXPECT_FALSE(b3->isNeighbor(b2));
    
    b1->removeNeighbor(b2);
    EXPECT_FALSE(b2->isNeighbor(b1));
    
    ///check compartment adding/removing
    Compartment* newC = GController::getCompartment(std::vector<double>{40,10,10});
    b1->setCompartment(newC);
    
    EXPECT_TRUE(newC->hasBoundaryElement(b1));
    EXPECT_FALSE(GController::getCompartment(std::vector<double>{15.0,15.0,15.0})->hasBoundaryElement(b1));
    
}

TEST(PlaneBoundaryElementTest, Distances) {
    
    SystemParameters::GParams.compartmentSizeX = 10.0;
    SystemParameters::GParams.compartmentSizeY = 10.0;
    SystemParameters::GParams.compartmentSizeZ = 10.0;
    
    SystemParameters::GParams.NX = 5;
    SystemParameters::GParams.NY = 5;
    SystemParameters::GParams.NZ = 5;
    
    SystemParameters::GParams.nDim = 3;
    
    GController::initializeGrid();
    
    BoundaryElement* b = new PlaneBoundaryElement({10.0,10.0,10.0}, {1,0,0}, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(5.0 ,b->distance({15,10,10}));
    EXPECT_EQ(-5.0 ,b->distance({5,10,10}));
    
    EXPECT_EQ(-5.0 ,b->stretchedDistance({15,15,10}, {1.0,0,0}, 10.0));
    EXPECT_EQ(-10.0 ,b->stretchedDistance({10,5,10}, {1.0,0,0}, 10.0));
}


#endif //DO_THIS_BOUNDARY_TEST