//
//  test_parser.cpp
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "gtest/gtest.h"
#include "Parser.h"

#define DO_THIS_PARSER_TEST

#ifdef DO_THIS_PARSER_TEST
#define TESTING


TEST(SystemParserTest, main) {
    
    SystemParser p("/Users/jameskomianos/Code/CytoSim-Repo/Cyto/TESTS/testsysteminput.txt");
    
    p.readBoundaryParameters();
    p.readGeometryParameters();
    p.readMechanicsParameters();
    
    EXPECT_TRUE(p.mechanics());
    EXPECT_TRUE(p.chemistry());
    
    ///Check a few
    EXPECT_EQ(10.0, SystemParameters::Boundaries().interactionCutoff);
    EXPECT_EQ(180.0, SystemParameters::Mechanics().FBendingTheta);
    EXPECT_EQ(0,SystemParameters::Mechanics().MStretchingL);
    EXPECT_EQ(100.0, SystemParameters::Geometry().compartmentSizeX);
    EXPECT_EQ(27.0, SystemParameters::Geometry().cylinderSize);
    
    ///Check string reading
    MechanicsAlgorithm MAlgorithm = p.readMechanicsAlgorithm();
    EXPECT_EQ(MAlgorithm.ConjugateGradient, "POLAKRIBIERE");
    
    ChemistryAlgorithm CAlgorithm = p.readChemistryAlgorithm();
    EXPECT_EQ(CAlgorithm.algorithm, "NRM");
    
    MechanicsFFType MType = p.readMechanicsFFType();
    EXPECT_EQ(MType.FBendingType, "HARMONIC");
    EXPECT_EQ(MType.FTwistingType, "");
    
    BoundaryType BType = p.readBoundaryType();
    EXPECT_EQ(BType.boundaryShape, "CUBIC" );
}

TEST(FilamentParserTest, main) {
    
    FilamentParser p("/Users/jameskomianos/Code/CytoSim-Repo/Cyto/TESTS/testfilamentinput.txt");
    
    auto filaments = p.readFilaments();
    std::vector<std::vector<std::vector<double>>> expectedValues = { {{10.0,10.0,10.0},{100.0,100.0,100.0}},
                                                                     {{50.0,50.0,50.0},{150.0,150.0,150.0}},
                                                                     {{100.0,100.0,100.0},{200.0,200.0,200.0}},
                                                                     {{300.0,300.0,300.0},{400.0,400.0,400.0}} };
    
    EXPECT_TRUE(std::equal(filaments.begin(), filaments.end(), expectedValues.begin()));
}




#endif ///DO_THIS_PARSER_TEST