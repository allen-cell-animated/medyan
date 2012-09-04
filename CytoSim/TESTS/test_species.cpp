//
//  test_species.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/17/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

// Note: This test omits many functions of Species that interact with Reaction objects. 
//        Separate tests weill cover those methods.

#define DO_THIS_TEST

#ifdef DO_THIS_TEST

#include <iostream>
#include "gtest/gtest.h"
#include "Species.h"
#include "SpeciesContainer.h"

using namespace std;
using namespace chem;

TEST(SpeciesNamesDBTest, All) {
    SpeciesNamesDB::Instance()->clear();
    int y = SpeciesNamesDB::Instance()->stringToInt("Arp2/3");
    EXPECT_EQ(0,y);
    std::string x = SpeciesNamesDB::Instance()->intToString(y);
    EXPECT_EQ("Arp2/3",x);
    
    EXPECT_THROW(SpeciesNamesDB::Instance()->intToString(1), std::out_of_range);

    y = SpeciesNamesDB::Instance()->stringToInt("G-Actin");
    EXPECT_EQ(1,y);
    x = SpeciesNamesDB::Instance()->intToString(y);
    EXPECT_EQ("G-Actin",x);
    EXPECT_NO_THROW(SpeciesNamesDB::Instance()->intToString(1));
}

TEST(SpeciesTest, CTors) {

    //Default CTor
    SpeciesBulk X;
    EXPECT_EQ(0,X.getN());
    EXPECT_EQ("",X.getName());
//    cout << X << endl;

    //CTor 1
    SpeciesBulk A{"Arp2/3",25};
    EXPECT_EQ(25,A.getN());
    EXPECT_EQ("Arp2/3",A.getName());
//    cout << A << endl;
    
    //Copy Ctor
    SpeciesBulk B(A);
    EXPECT_EQ(A,B);
//    cout << B << endl;

    
    //Move Assignment
    SpeciesBulk D{"G-Actin",300};
//    cout << D << endl;
    D = [] () {SpeciesBulk Y("G-Actin",39); return Y;}();
//    cout << D << endl;
    EXPECT_EQ(39,D.getN());
    EXPECT_EQ("G-Actin",D.getName());
    
    //Assignment operator
    SpeciesBulk F{"Myosin-X",112};
//    cout << F << endl;
    F=A;
    EXPECT_EQ(A,F);

    
}

TEST(SpeciesTest, CopyNumber) {
    SpeciesBulk B{"Arp2/3",25};
    B.setN(32);
    EXPECT_EQ(32,B.getN());
}


TEST(SpeciesTest, Iostream) {    
    SpeciesBulk B{"Arp2/3",25};

    ostringstream outStream;
    
    outStream << B;
    string res = outStream.str();
    string b_fullname = B.getFullName() + "[" + std::to_string(B.getN()) + "]";
    EXPECT_EQ(b_fullname,res);
}

TEST(SpeciesTest, Vector) {
    vector<SpeciesBulk> vsb;
    vsb.push_back({"G-Actin",133});
//    cout << vsb.back() << endl;
    EXPECT_EQ(133,vsb.back().getN());
    EXPECT_EQ("G-Actin",vsb.back().getName());
    
    vsb.emplace_back("Arp2/3",244);
//    cout << vsb.back() << endl;
    EXPECT_EQ("Arp2/3",vsb.back().getName());

    EXPECT_NE(vsb[0],vsb[1]);
    
    vsb.push_back({"Motor",800});
    EXPECT_EQ("Motor",vsb[2].getName());
    RSpecies &rs_before = vsb[2].getRSpecies();
    
    for(int i=0; i<10; ++i){
//        cout << i << ":" << endl;
        vsb.push_back({"Motor",900});
//        cout << vsb.back() << endl;
    }
    
    /// To make sure that move should be enabled:
    EXPECT_TRUE(is_nothrow_move_constructible<SpeciesBulk>::value);
    EXPECT_TRUE(is_copy_constructible<SpeciesBulk>::value);
    
    /// Now checking if Species have been moved and not copied by the vector<Species> reallocators
    RSpecies &rs_after = vsb[2].getRSpecies();
    EXPECT_EQ(&rs_before,&rs_after); // This makes sure that the move constructors were involved during the internal 
                                     // reallocation of the vector<Species>, and RSpecies pointer was conserved
    
//    cout << "Type Traits: " << boolalpha << is_nothrow_move_constructible<SpeciesBulk>::value << ", " 
//    << is_copy_constructible<SpeciesBulk>::value << endl;
//
//    cout << "Still in TEST(SpeciesTest, Vector) {..." << endl;
}

TEST(SpeciesContainerTest, SpeciesPtrContainerVector) {
    SpeciesPtrContainerVector scv;
    
    Species *x0 = scv.addSpecies<SpeciesBulk>("Profilin",31);
    Species *x1 = scv.addSpecies<SpeciesDiffusing>("Actin",99);
    Species *x2 = scv.addSpecies<SpeciesDiffusing>("Arp2/3",11);
    Species *x3 = scv.addSpecies<SpeciesDiffusing>("Capping",22);
    
    EXPECT_EQ(x0, scv.findSpeciesByIndex(0));
    EXPECT_EQ(x1, scv.findSpeciesByIndex(1));
    EXPECT_EQ(x2, scv.findSpeciesByIndex(2));
    EXPECT_EQ(x3, scv.findSpeciesByIndex(3));

    EXPECT_EQ("Capping", scv.findSpeciesByIndex(3)->getName());

    //cout << (*x1) << (*x2) << endl;
    //scv.printSpecies();
    //cout << endl;
    
    scv.removeSpecies(x0);
    EXPECT_EQ(x3, scv.findSpeciesByIndex(2));
//    scv.printSpecies();
//    cout << endl;
    
    scv.removeSpecies("Actin");
//    scv.printSpecies();
    EXPECT_EQ(x3, scv.findSpeciesByIndex(1));
//    cout << endl;
    
    Species *y = scv.findSpeciesByName("Arp2/3");
    EXPECT_EQ(x2, y);
//    cout << (*y) << endl;
}

TEST(SpeciesContainerTest, SpeciesContainerVector) {
    SpeciesContainerVector<SpeciesDiffusing> scv;
    
    size_t profilin = scv.addSpecies("Profilin",species_copy_t(31));
    size_t actin = scv.addSpecies("Actin",species_copy_t(99));
    size_t arp23 = scv.addSpecies("Arp2/3",species_copy_t(11));
    size_t capping = scv.addSpecies("Capping",species_copy_t(22));
    
    EXPECT_EQ(size_t(0), profilin);
    EXPECT_EQ(size_t(1), actin);
    EXPECT_EQ(size_t(2), arp23);
    EXPECT_EQ(size_t(3), capping);
    
    EXPECT_EQ("Profilin", scv.findSpecies(0).getName());
    EXPECT_EQ("Actin", scv.findSpecies(1).getName());
    EXPECT_EQ("Arp2/3", scv.findSpecies(2).getName());
    EXPECT_EQ("Capping", scv.findSpecies(3).getName());
    
    scv.removeSpecies("Profilin");
    EXPECT_EQ("Capping", scv.findSpecies(2).getName());
    //    scv.printSpecies();
    
    scv.removeSpecies("Actin");
    //    scv.printSpecies();
    EXPECT_EQ("Capping", scv.findSpecies(1).getName());
    
    Species &y = scv.findSpecies("Arp2/3");
    EXPECT_EQ("Arp2/3", y.getName());
    
    size_t index = scv.findSpeciesIndex("Arp2/3");
    EXPECT_EQ(0,index);
}

#endif // DO_THIS_TEST