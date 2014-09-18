//
//  test_geometry.cpp
//  Cyto
//
//  Created by James Komianos on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

//#define DO_THIS_GEOMETRY_TEST

#ifdef DO_THIS_GEOMETRY_TEST

#include <iostream>
#include "gtest/gtest.h"

#include "common.h"
#include "GController.h"

using namespace std;

///testing basic initialization and getters
TEST(GeometryTest, Main) {

    SystemParameters::GParams.compartmentSizeX = 10.0;
    SystemParameters::GParams.compartmentSizeY = 10.0;
    SystemParameters::GParams.compartmentSizeZ = 10.0;
    
    SystemParameters::GParams.NX = 5;
    SystemParameters::GParams.NY = 5;
    SystemParameters::GParams.NZ = 5;
    
    SystemParameters::GParams.monomerSize = 2.7;
    SystemParameters::GParams.cylinderSize = 27.0;
    
    SystemParameters::GParams.nDim = 3;
    
    GController::initializeGrid();
    
    ///checking out of bounds
    EXPECT_ANY_THROW(GController::getCompartment(std::vector<double>{60.0,50.0,50.0}));
    
    EXPECT_EQ(GController::getCompartment(std::vector<size_t>{0,0,0}), GController::getCompartment(std::vector<double>{5.0,5.0,5.0}));
    EXPECT_EQ(GController::getCompartment(std::vector<size_t>{0,1,0}), GController::getCompartment(std::vector<double>{5.0,15.0,5.0}));
    EXPECT_EQ(GController::getCompartment(std::vector<size_t>{0,1,1}), GController::getCompartment(std::vector<double>{5.0,15.0,15.0}));
}


#endif //DO_THIS_GEOMETRY_TEST
