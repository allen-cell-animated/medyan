//
//  test_parser.cpp
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

//#define DO_THIS_PARSER_TEST

#ifdef DO_THIS_PARSER_TEST

#include "gtest/gtest.h"

#include "common.h"
#include "Parser.h"
#include "SystemParameters.h"

TEST(SystemParserTest, Main) {
    
    SystemParser p("/Users/jimmy/Code/Cyto/Cyto/testsysteminput.txt");
    
    p.readBoundaryParameters();
    p.readGeometryParameters();
    p.readMechanicsParameters();
    
    ///Check a few
    EXPECT_EQ(10.0, SystemParameters::Boundaries().boundaryCutoff);
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

TEST(FilamentParserTest, Main) {
    
    FilamentParser p("/Users/jimmy/Code/Cyto/Cyto/testfilamentinput.txt");
    
    auto filaments = p.readFilaments();
    vector<vector<vector<double>>> expectedValues = { {{100.0,100.0,100.0},{200.0,200.0,200.0}}};
    
    EXPECT_TRUE(equal(filaments.begin(), filaments.end(), expectedValues.begin()));
}

TEST(ChemParserTest, Main) {
    
    ChemistryParser p("/Users/Konstantin/Documents/Codes/Cyto/CytoRepo/Cyto/testchemistryinput.txt");
    
    ChemistryData chem = p.readChemistryInput();
    EXPECT_EQ("A", get<0>(chem.speciesBulk[0]));
    EXPECT_EQ(1000, get<1>(chem.speciesBulk[0]));
    
    EXPECT_EQ("A", chem.speciesFilament[0]);
    EXPECT_EQ("E", chem.speciesBound[0]);
    EXPECT_EQ("P", chem.speciesPlusEnd[0]);
    EXPECT_EQ("M", chem.speciesMinusEnd[0]);
    
    EXPECT_EQ("A:BULK", get<0>(chem.reactions[0])[0]);
    EXPECT_EQ("P:PLUSEND:N", get<0>(chem.reactions[0])[1]);
    EXPECT_EQ("A:FILAMENT", get<1>(chem.reactions[0])[0]);
    EXPECT_EQ("E:BOUND", get<1>(chem.reactions[0])[1]);
    EXPECT_EQ("P:PLUSEND:N+1", get<1>(chem.reactions[0])[2]);
    EXPECT_EQ(10.0, get<2>(chem.reactions[0]));
    
    EXPECT_EQ("A:BULK", get<0>(chem.reactions[1])[0]);
    EXPECT_EQ("M:MINUSEND:N+1", get<0>(chem.reactions[1])[1]);
    EXPECT_EQ("A:FILAMENT", get<1>(chem.reactions[1])[0]);
    EXPECT_EQ("E:BOUND", get<1>(chem.reactions[1])[1]);
    EXPECT_EQ("M:MINUSEND:N", get<1>(chem.reactions[1])[2]);
    EXPECT_EQ(10.0, get<2>(chem.reactions[1]));
    
    EXPECT_EQ("A:BULK", get<1>(chem.reactions[2])[0]);
    EXPECT_EQ("P:PLUSEND:N", get<1>(chem.reactions[2])[1]);
    EXPECT_EQ("A:FILAMENT", get<0>(chem.reactions[2])[0]);
    EXPECT_EQ("E:BOUND", get<0>(chem.reactions[2])[1]);
    EXPECT_EQ("P:PLUSEND:N+1", get<0>(chem.reactions[2])[2]);
    EXPECT_EQ(10.0, get<2>(chem.reactions[2]));
    
}




#endif ///DO_THIS_PARSER_TEST