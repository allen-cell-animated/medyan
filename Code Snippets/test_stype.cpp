//
//  test_stype.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/17/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "gtest/gtest.h"

#include "SpeciesType.h"

using namespace std;
using namespace boost::flyweights;
using namespace chem;


TEST(SpeciesTypeTest, CTors) {
    // Testing the main ctor
    SpeciesType A("Arp2/3",SType::Diffusing);
    EXPECT_EQ("Arp2/3", A.getName());
    EXPECT_EQ(SType::Diffusing, A.getType());
    
    // Testing the copy ctor
    SpeciesType B(A);
    EXPECT_EQ(A,B);
}

TEST(SpeciesTypeTest, Assignments) {
    SpeciesType A("Arp2/3",SType::Diffusing);
    SpeciesType C("X",SType::Bulk);
    EXPECT_EQ("X", C.getName());
    C=A;
    EXPECT_EQ(A,C);
    
    SpeciesType D("D",SType::Membrane);
    EXPECT_NE(C,D);
}

TEST(SpeciesTypeTest, Accessors) {
    SpeciesType A("Arp2/3",SType::Diffusing);
    SpeciesType C("X",SType::Bulk);

    EXPECT_EQ("Arp2/3", A.getName());
    EXPECT_EQ("X", C.getName());
    
    EXPECT_EQ("Diffusing",A.getTypeAsString());
    EXPECT_EQ("Bulk",C.getTypeAsString());
    
    EXPECT_EQ(SType::Diffusing, A.getType());
    EXPECT_EQ(SType::Bulk, C.getType());
    
    EXPECT_TRUE(A.is_of_type("Arp2/3",SType::Diffusing));
    EXPECT_TRUE(C.is_of_type("X",SType::Bulk));
}

TEST(SpeciesTypeTest, InsertingIntoVector) {
    vector<SpeciesType> vst;

    vst.push_back({"Arp2/3",SType::Diffusing});
    EXPECT_EQ("Arp2/3",vst[0].getName());

    for(int i=0; i<100000; ++i){
        vst.emplace_back("X",SType::Bulk);
    }
    
    EXPECT_EQ("X",vst.back().getName());
    EXPECT_EQ("Bulk",vst.back().getTypeAsString());
    EXPECT_EQ(SType::Bulk,vst.back().getType());
    
    auto num_elem = vst.size();
    EXPECT_EQ(vst[num_elem-1],vst[num_elem-2]);
}

TEST(SpeciesTypeTest, ReturningFromAFunction) {
    auto f = [](){return SpeciesType("X",SType::Bulk);}; // A C++11 lambda f-n
    auto Y = f(); // calling it
    
    SpeciesType Z("X",SType::Bulk);
    EXPECT_EQ(Y,Z);
}

TEST(SpeciesTypeTest, Iostream) {
    SpeciesType B("Arp2/3",SType::Diffusing);
    ostringstream outStream;
    
    outStream << B;
    string res = outStream.str();
    string b_fullname = B.getName() + "{" + B.getTypeAsString() + "}";
    EXPECT_EQ(b_fullname,res);
}

TEST(SpeciesTypeTest, BoostFlywieght) {
    SpeciesType A = {"Arp2/3",SType::Diffusing};
    flyweight<SpeciesType> fA(A);

    flyweight<SpeciesType> fB (SpeciesType("Arp2/3",SType::Diffusing));
    EXPECT_EQ(fA.get(), fB.get());
    
    SpeciesType B(fB); // testing conversion from flyweight<SpeciesType> to SpeciesType
    EXPECT_EQ(B,A);
    
    ostringstream out1;
    out1 << B;
    ostringstream out2;
    out2 << fB;
    EXPECT_EQ(out1.str(), out2.str());
    
    // check hashing
    auto a = hash_value(A);
    auto b = hash_value(B);
    EXPECT_EQ(a,b);
    SpeciesType C = {"Arp2/3",SType::Bulk};
    auto c = hash_value(C);
    EXPECT_NE(a,c);
    
    vector<flyweight<SpeciesType>> vfst;
    for(int i=0; i<100; ++i){
        vfst.emplace_back("X",SType::Bulk);
    }
    const SpeciesType &m = vfst[10].get();
    const SpeciesType &n = vfst[90].get();
    EXPECT_EQ(&m,&n);
}

TEST(SpeciesTypeTest, BoostSerializationSerial) {
    
    SpeciesType *AA1 = new SpeciesType{"Arp2/3",SType::Bulk};
    SpeciesType *AA2 = nullptr;

    // create and open a character archive for output
    {
        ofstream ofs{"filename_x1y1z1.txt",ios::out};
        boost::archive::text_oarchive oa(ofs); 
        oa << AA1;
    }
    
    // read back the object
    {
        ifstream ifs{"filename_x1y1z1.txt",ios::in};
        boost::archive::text_iarchive ia(ifs);
        ia >> AA2;
//        cout << "BoostSerialization: {" << *AA2 << "}" << endl;
    }
    EXPECT_EQ(*AA1,*AA2);
}

TEST(SpeciesTypeTest, BoostSerializationVector) {
    std::vector<SpeciesType> vst;    
    for(int i=0; i<100; ++i){
        vst.emplace_back("X",SType::Bulk);
    }
    // create and open a character archive for output
    {
        ofstream ofs{"filename_x2y2z2.txt",ios::out};
        boost::archive::text_oarchive oa(ofs); 
        oa & vst;
    }
    
    // read back the object
    std::vector<SpeciesType> vst2;    
    {
        ifstream ifs{"filename_x2y2z2.txt",ios::in};
        boost::archive::text_iarchive ia(ifs);
        ia >> vst2;
        //        cout << "BoostSerialization: {" << *AA2 << "}" << endl;
    }
    
    EXPECT_EQ(vst[0],vst2[0]);    
    EXPECT_EQ(vst.back(),vst2.back());    
}

