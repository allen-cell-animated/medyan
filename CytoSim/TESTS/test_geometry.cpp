//
//  test_geometry.cpp
//  Cyto
//
//  Created by James Komianos on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#define DO_THIS_GEOMETRY_TEST

#ifdef DO_THIS_GEOMETRY_TEST
#define TESTING

#include <iostream>
#include "gtest/gtest.h"
#include "GController.h"
#include "CompartmentContainer.h"

using namespace std;

///testing basic initialization
TEST(GeometryTest, Initialization) {
    
    
    GController g;
    g.initializeGrid(3, {3,3,3}, {100.0,100.0,100.0});
    
    EXPECT_EQ(GController::_nDim, 3);
    
    
    
    EXPECT_EQ(GController::_systemSize[0], 300.0);
    EXPECT_EQ(GController::_systemSize[0], 300.0);
    EXPECT_EQ(GController::_systemSize[0], 300.0);
    
    EXPECT_EQ(GController::_compartmentSize[0], 100.0);
}








#endif //DO_THIS_GEOMETRY_TEST
