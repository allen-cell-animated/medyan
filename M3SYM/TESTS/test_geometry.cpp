
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifdef TESTING

//#define DO_THIS_GEOMETRY_TEST
#ifdef DO_THIS_GEOMETRY_TEST

#include "gtest/gtest.h"

#include "common.h"

#include "GController.h"
#include "SysParams.h"

///testing basic initialization and getters
TEST(GeometryTest, Basic) {

    SysParams::GParams.compartmentSizeX = 10.0;
    SysParams::GParams.compartmentSizeY = 10.0;
    SysParams::GParams.compartmentSizeZ = 10.0;
    
    SysParams::GParams.NX = 5;
    SysParams::GParams.NY = 5;
    SysParams::GParams.NZ = 5;
    
    SysParams::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    ///checking out of bounds
    EXPECT_ANY_THROW(GController::getCompartment(vector<double>{60.0,50.0,50.0}));
    
    //checking equivalence
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,0,0}),
              GController::getCompartment(vector<double>{5.0,5.0,5.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,1,0}),
              GController::getCompartment(vector<double>{5.0,15.0,5.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,1,1}),
              GController::getCompartment(vector<double>{5.0,15.0,15.0}));
    
    //checking edges
    EXPECT_EQ(GController::getCompartment(vector<size_t>{4,4,4}),
              GController::getCompartment(vector<double>{49.0,49.0,49.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{4,0,0}),
              GController::getCompartment(vector<double>{49.0,0.0,0.0}));
    EXPECT_EQ(GController::getCompartment(vector<size_t>{0,0,0}),
              GController::getCompartment(vector<double>{0.0,0.0,0.0}));
    
}

TEST(GeometryTest, NonCubicGrid) {
    
    SysParams::GParams.compartmentSizeX = 10.0;
    SysParams::GParams.compartmentSizeY = 20.0;
    SysParams::GParams.compartmentSizeZ = 100.0;
    
    SysParams::GParams.NX = 15;
    SysParams::GParams.NY = 5;
    SysParams::GParams.NZ = 10;
    
    SysParams::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
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
#endif //TESTING
