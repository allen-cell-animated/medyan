
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifdef TESTING

#define DO_THIS_GEOMETRY_TEST
#ifdef DO_THIS_GEOMETRY_TEST

#include "gtest/gtest.h"

#include "common.h"

#include "GController.h"
#include "SysParams.h"
///testing basic initialization and getters
TEST(GeometryTest, Basic) {

    SysParams
//  http://papoian.chem.umd.edu/::GParams.compartmentSizeX = 10.0;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.compartmentSizeY = 10.0;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.compartmentSizeZ = 10.0;
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.NX = 5;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.NY = 5;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.NZ = 5;
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.monomerSize = 2.7;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.cylinderSize = 27.0;
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
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
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.compartmentSizeX = 10.0;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.compartmentSizeY = 20.0;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.compartmentSizeZ = 100.0;
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.NX = 15;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.NY = 5;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.NZ = 10;
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.monomerSize = 1;
    SysParams
//  http://papoian.chem.umd.edu/::GParams.cylinderSize = 5.0;
    
    SysParams
//  http://papoian.chem.umd.edu/::GParams.nDim = 3;
    
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
