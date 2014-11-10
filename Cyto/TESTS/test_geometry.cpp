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
#include "SystemParameters.h"

///testing basic initialization and getters
TEST(GeometryTest, Basic) {

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
    EXPECT_ANY_THROW(GController::getCompartment(vector<double>{60.0,50.0,50.0}));
    
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,0,0}),
              GController::getCompartment(vector<double>{5.0,5.0,5.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,1,0}),
              GController::getCompartment(vector<double>{5.0,15.0,5.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,1,1}),
              GController::getCompartment(vector<double>{5.0,15.0,15.0}));
}

TEST(GeometryTest, NonCubicGrid) {
    
    SystemParameters::GParams.compartmentSizeX = 10.0;
    SystemParameters::GParams.compartmentSizeY = 20.0;
    SystemParameters::GParams.compartmentSizeZ = 100.0;
    
    SystemParameters::GParams.NX = 15;
    SystemParameters::GParams.NY = 5;
    SystemParameters::GParams.NZ = 10;
    
    SystemParameters::GParams.monomerSize = 1;
    SystemParameters::GParams.cylinderSize = 5.0;
    
    SystemParameters::GParams.nDim = 3;
    
    GController::initializeGrid();
    
    EXPECT_ANY_THROW(GController::getCompartment(vector<double>{60.0,50.0,1050.0}));
    EXPECT_ANY_THROW(GController::getCompartment(vector<double>{200.0,50.0,900.0}));
    EXPECT_ANY_THROW(GController::getCompartment(vector<double>{100.0,110.0,900.0}));
    
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,0,0}),
              GController::getCompartment(vector<double>{5.0,5.0,5.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,1,0}),
              GController::getCompartment(vector<double>{5.0,25.0,5.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,1,1}),
              GController::getCompartment(vector<double>{5.0,30.0,190.0}));
    
}



#endif //DO_THIS_GEOMETRY_TEST
