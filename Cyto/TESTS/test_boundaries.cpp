//
//  test_boundaries.cpp
//  Cyto
//
//  Created by James Komianos on 9/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#define DO_THIS_BOUNDARY_TEST

#ifdef DO_THIS_BOUNDARY_TEST

#include <iostream>
#include "gtest/gtest.h"

#include "common.h"
#include "BoundaryImpl.h"
#include "SystemParameters.h"

TEST(PlaneBoundaryElementTest, Distances) {
    
    SystemParameters::GParams.compartmentSizeX = 10.0;
    SystemParameters::GParams.compartmentSizeY = 10.0;
    SystemParameters::GParams.compartmentSizeZ = 10.0;
    
    SystemParameters::GParams.NX = 5;
    SystemParameters::GParams.NY = 5;
    SystemParameters::GParams.NZ = 5;
    
    SystemParameters::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    BoundaryElement* b = new PlaneBoundaryElement({10.0,10.0,10.0}, {1,0,0}, 1.0, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(5.0 ,b->distance({15,10,10}));
    EXPECT_EQ(-5.0 ,b->distance({5,10,10}));
    
    EXPECT_EQ(15.0 ,b->stretchedDistance({15,15,10}, {1.0,0,0}, 10.0));
    EXPECT_EQ(10.0 ,b->stretchedDistance({10,5,10}, {1.0,0,0}, 10.0));
}

TEST(SphereBoundaryElementTest, Distances) {
    
    SystemParameters::GParams.compartmentSizeX = 10.0;
    SystemParameters::GParams.compartmentSizeY = 10.0;
    SystemParameters::GParams.compartmentSizeZ = 10.0;
    
    SystemParameters::GParams.NX = 5;
    SystemParameters::GParams.NY = 5;
    SystemParameters::GParams.NZ = 5;
    
    SystemParameters::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    BoundaryElement* b = new SphereBoundaryElement({25.0,25.0,25.0}, 10.0, 1.0, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(10.0, b->distance({25.0,25.0,25.0}));
    EXPECT_EQ(5.0, b->distance({25.0,30.0,25.0}));
    EXPECT_EQ(-5.0, b->distance({25.0, 40.0, 25.0}));
    
    EXPECT_EQ(0 ,b->stretchedDistance({25,25,25}, {1.0,0,0}, 10.0));
    EXPECT_EQ(-10.0 ,b->stretchedDistance({25,35,25}, {0,1,0}, 10.0));

}

TEST(CubicBoundary, Within) {
    
    
    
    
    
}






#endif //DO_THIS_BOUNDARY_TEST