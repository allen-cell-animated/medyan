//
//  test_simple_gillespie.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 8/18/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

//#define DO_THIS_TEST

#ifdef DO_THIS_TEST


#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>

using namespace boost::accumulators;

#include "gtest/gtest.h"

#include "Species.h"
#include "Reaction.h"
#include "ChemSimpleGillespieImpl.h"
#include "ChemSim.h"

using namespace std;
using namespace chem;

TEST(ChemSimpleGillespieTest, SimpleStoichiometryInvariants) {
    SpeciesBulk A1("A1",  100);
    SpeciesBulk A2("A2", 0);
    Reaction<1,1> r1 = { {&A1,&A2}, 100.0 };
    
    ChemSimpleGillespieImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.initialize();
    //    chem.printReactions();
    chem.run(30);
    //    chem.printReactions();
    EXPECT_EQ(70,A1.getN());
    EXPECT_EQ(30,A2.getN());
}

TEST(ChemSimpleGillespieTest, StoichiometryInvariants) {
    SpeciesBulk A1("A1",  100);
    SpeciesBulk A2("A2", 0);
    SpeciesBulk A3("A3", 0);
    Reaction<1,1> r1 = { {&A1,&A2}, 10.0 };
    Reaction<1,1> r2 = { {&A2,&A1}, 15.0 };
    Reaction<1,1> r3 = { {&A1,&A3}, 20.0 };
    
    ChemSimpleGillespieImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.addReaction(&r3);
    chem.initialize();
    //    chem.printReactions();
    chem.run(30);
    //    chem.printReactions();
    EXPECT_EQ(100,A1.getN()+A2.getN()+A3.getN());
}

// Testing A<->B steady state
TEST(ChemSimpleGillespieTest, SimpleSteadyState) {
    int Nstart = 16;
    SpeciesBulk A1("A1",  Nstart);
    SpeciesBulk A2("A2", 0);
    // A1 <-> A2 with the same forward and backward rates; [A]~[B] at steady state
    Reaction<1,1> r1 = { {&A1,&A2}, 100.0 };
    Reaction<1,1> r2 = { {&A2,&A1}, 100.0 };
    
    ChemSimpleGillespieImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.initialize();
    //    chem.printReactions();
    chem.run(10000);
    // chem.printReactions();
    accumulator_set<int, stats<tag::variance(immediate)>> accA1;
    accumulator_set<int, stats<tag::mean>> accA2;
    accumulator_set<int, stats<tag::covariance<double, tag::covariate1> > > accCov;
    int N_SAMPLE_POINTS=1000;
    for(int i=0;i<N_SAMPLE_POINTS;++i){
        chem.run(100);
        //        cout << "N= " << A1.getN() << ", " << A1.getN() << endl;
        accA1(A1.getN());
        accCov(A1.getN(), covariate1 = A2.getN());
    }
    double A1mean = mean(accA1);
    double A1var = variance(accA1);
    double var_expected = double(Nstart)/2/2;
    double mean_error = sqrt(var_expected/N_SAMPLE_POINTS);
    EXPECT_TRUE(abs(A1mean-Nstart/2)<10*mean_error);
    EXPECT_TRUE(abs(var_expected-A1var)<0.1*var_expected);// within 10% of the expected variance
    EXPECT_FLOAT_EQ(-1.0,covariance(accCov)/(A1var));
    //    cout << "ChemSimpleGillespieTest::SimpleSteadyState: " << A1mean << " vs " << double(Nstart)/2 << ", Expected error: " << mean_error << endl;
    //    cout << "ChemSimpleGillespieTest::SimpleSteadyState: std= " << sqrt(A1var) << ", "  << covariance(accCov)/A1var << endl;
}

// Testing transient dynamics for A<->B
TEST(ChemSimpleGillespieTest, SimpleTransient) {
    const long long int N_SAMPLE_POINTS=pow(10,6);
    const long long int Nstart = 10;
    const double tau_snapshot = 0.48; //seconds
    //long long int print_freq = pow(10,7);
    SpeciesBulk A1("A1",  Nstart);
    SpeciesBulk A2("A2", 0);
    // A1 <-> A2 with the same forward and backward rates; [A]~[B] at steady state
    Reaction<1,1> r1 = { {&A1,&A2}, 2.5 };
    Reaction<1,1> r2 = { {&A2,&A1}, 2.5 };
    
    ChemSimpleGillespieImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    //    chem.printReactions();
    // chem.printReactions();
    vector<long long int> n_hist(Nstart+1);
    
    //int counter = 0;
    accumulator_set<double, stats<tag::mean>> accTau;
    long long int N_penultimate;
    for(long long int i=0;i<N_SAMPLE_POINTS;++i){
        A1.setN(Nstart);
        A2.setN(0);
        chem.initialize();
        //        long long int events=0;
        do {
            N_penultimate=A1.getN();
            chem.run(1);
            //            ++events;
        } while (tau()<tau_snapshot);
        ++n_hist[N_penultimate];
        //        ++n_hist[A1.getN()];
        accTau(tau());
        
        //        if(i%print_freq==0)
        //            cout << "i=" << i << endl;
    }
    
    //    cout << "tau_mean=" << mean(accTau) << ", counter=" << N_SAMPLE_POINTS << endl;;
    //    for (long long int n=0; n<=Nstart; ++n){
    //        cout << "P[" << n << "]=" << double(n_hist[n])/N_SAMPLE_POINTS << endl;
    //    }
    
    double sum=0;
    for(auto num: n_hist)
        sum+=double(num)/N_SAMPLE_POINTS;
    
    // The results below are for N=10 (coming from both analytical formula and numerical integration)
    vector<double> n_hist_analyt {0.0003773 ,  0.00452585,  0.02443017,  0.0781464 ,  0.1640442 , 0.23613261,  0.23604161,  0.1617947,  0.07277957,  0.01940041, 0.00232715};
    double relative_error=0.15; //i.e. allow a 15% relative error
    for(int n=0; n<(Nstart+1); ++n){
        double p_est=double(n_hist[n])/N_SAMPLE_POINTS;
        double p_analyt=n_hist_analyt[n];
        //        cout << "P[" << n << "]=" << p_est << " " << p_analyt << endl;
        EXPECT_NEAR(p_est,p_analyt,relative_error*p_analyt);
    }
    // Note that the error for P[0] is hard to converge, most likely because of P[0] being very small
}


// Testing transient dynamics for the A->B->C cycle, where A, B, and C can only take two values, n=0,1
TEST(ChemSimpleGillespieTest, CyclicTransient) {
    const long long int N_SAMPLE_POINTS=pow(10,6);
    const double tau_snapshot = 0.25; //seconds
    //long long int print_freq = pow(10,7);
    SpeciesBulk A1("A1", 1, 1);
    SpeciesBulk A2("A2", 0, 1);
    SpeciesBulk A3("A3", 0, 1);
    
    // cout << "A2: upper limit: " << A2.getUpperLimitForN() << endl;
    // A1 <-> A2 with the same forward and backward rates; [A]~[B] at steady state
    Reaction<1,1> r1 = { {&A1,&A2}, 4.5 };
    Reaction<1,1> r2 = { {&A2,&A3}, 2.5 };
    Reaction<1,1> r3 = { {&A3,&A1}, 0.5 };
    
    ChemSimpleGillespieImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.addReaction(&r3);
    
    long long int n_a1_hist=0;
    long long int n_a2_hist=0;
    long long int n_a3_hist=0;
    
    for(long long int i=0;i<N_SAMPLE_POINTS;++i){
        A1.setN(1);
        A2.setN(0);
        A3.setN(0);
        long long int n_a1_pentult=0;
        long long int n_a2_pentult=0;
        long long int n_a3_pentult=0;
        chem.initialize();
        long long int events=0;
        do {
            n_a1_pentult=A1.getN();
            n_a2_pentult=A2.getN();
            n_a3_pentult=A3.getN();
            chem.run(1);
            //            cout << "Event [" << events << "], tau=" << tau() << ",  A1=" << A1.getN() << ", A2= " << A2.getN() << ", A3= " << A3.getN() << endl;
            ++events;
            //            chem.printReactions();
            //            cout << endl << endl;
        } while (tau()<tau_snapshot);
        n_a1_hist+=n_a1_pentult;
        n_a2_hist+=n_a2_pentult;
        n_a3_hist+=n_a3_pentult;
        //        cout << "1=" << n_a1_hist << ", 2=" << n_a2_hist << ", 3=" << n_a3_hist << endl;
    }
    
    double pa1=static_cast<double>(n_a1_hist)/N_SAMPLE_POINTS;
    double pa2=static_cast<double>(n_a2_hist)/N_SAMPLE_POINTS;
    double pa3=static_cast<double>(n_a3_hist)/N_SAMPLE_POINTS;
    
    double pa1_numeric = 0.33169986;
    double pa2_numeric = 0.47589009;
    double pa3_numeric = 0.19241006;
    
    //    cout << "P1=" << pa1 << ", P1_numeric=" << pa1_numeric << endl
    //        << "P2=" << pa2 << ", P2_numeric=" << pa2_numeric << endl
    //        << "P3=" << pa3 << ", P3_numeric=" << pa3_numeric << endl;
    
    double relative_error=0.01; //i.e. allow a 1% relative error
    EXPECT_NEAR(pa1,pa1_numeric,relative_error*pa1_numeric);
    EXPECT_NEAR(pa2,pa2_numeric,relative_error*pa1_numeric);
    EXPECT_NEAR(pa3,pa3_numeric,relative_error*pa1_numeric);
    
}

// Testing transient dynamics for the X<->A->B->C cycle, where A, B, and C can only take two values, n=0,1
#ifdef TRACK_UPPER_COPY_N
TEST(ChemSimpleGillespieTest, ComplexCyclicTransient) {
    const long long int N_SAMPLE_POINTS=pow(10,6);
    const long long int Nstart = 3;
    const double tau_snapshot = 0.5; //seconds
    //long long int print_freq = pow(10,7);
    SpeciesBulk X("X", Nstart); // X's copy number is not restricted
    SpeciesBulk A("A", 0, 1);
    SpeciesBulk B("B", 0, 1);
    SpeciesBulk C("C", 0, 1);
    
    // cout << "A2: upper limit: " << A2.getUpperLimitForN() << endl;
    // A1 <-> A2 with the same forward and backward rates; [A]~[B] at steady state
    float kxa=0.8; // s^-1
    float kax=0.3; // s^-1
    float kab=4.5; // s^-1
    float kbc=2.5; // s^-1
    float kca=0.5; // s^-1
    Reaction<1,1> xa = { {&X,&A}, kxa };
    Reaction<1,1> ax = { {&A,&X}, kax };
    Reaction<1,1> r1 = { {&A,&B}, kab };
    Reaction<1,1> r2 = { {&B,&C}, kbc };
    Reaction<1,1> r3 = { {&C,&A}, kca };
    
    ChemSimpleGillespieImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&xa);
    chem.addReaction(&ax);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.addReaction(&r3);
    
    vector<long long int> x_hist(Nstart+1);
    long long int n_a1_hist=0;
    long long int n_a2_hist=0;
    long long int n_a3_hist=0;
    
    for(long long int i=0;i<N_SAMPLE_POINTS;++i){
        A.setN(0);
        B.setN(0);
        C.setN(0);
        X.setN(Nstart);
        long long int n_a_pentult=0;
        long long int n_b_pentult=0;
        long long int n_c_pentult=0;
        long long int x_pentult=0;
        chem.initialize();
        long long int events=0;
        do {
            //            cout << "chem.run(1) start, i= " << i << endl;
            x_pentult=X.getN();
            n_a_pentult=A.getN();
            n_b_pentult=B.getN();
            n_c_pentult=C.getN();
            bool success = chem.run(1);
            if(!success){
                cout << "chem.run(1) has failed, i= " << i << endl;
                chem.printReactions();
                break;
            }
            //            cout << "Event [" << events << "], tau=" << tau() << ",  A1=" << A1.getN() << ", A2= " << A2.getN() << ", A3= " << A3.getN() << endl;
            ++events;
            //            chem.printReactions();
            //            cout << endl << endl;
        } while (tau()<tau_snapshot);
        ++x_hist[x_pentult];
        n_a1_hist+=n_a_pentult;
        n_a2_hist+=n_b_pentult;
        n_a3_hist+=n_c_pentult;
        //        cout << "1=" << n_a1_hist << ", 2=" << n_a2_hist << ", 3=" << n_a3_hist << endl;
    }
    
    
    vector<double> p_nrm;
    
    for(int n=0; n<(Nstart+1); ++n){
        double p_est=double(x_hist[n])/N_SAMPLE_POINTS;
        p_nrm.push_back(p_est);
    }
    
    double pa1=static_cast<double>(n_a1_hist)/N_SAMPLE_POINTS;
    double pa2=static_cast<double>(n_a2_hist)/N_SAMPLE_POINTS;
    double pa3=static_cast<double>(n_a3_hist)/N_SAMPLE_POINTS;
    p_nrm.push_back(pa1);
    p_nrm.push_back(pa2);
    p_nrm.push_back(pa3);
    
    // The results below are for ...
    vector<double> p_numeric {0.001687323512088279, 0.12264078507458409, 0.55515007879166167, 0.3205218126216664, 0.32672439967797662, 0.30766594955383336, 0.17110327024528463};
    double relative_error=0.05; //i.e. allow a 5% relative error
    
    for(int n=0; n<(Nstart+4); ++n){
        //        cout << "P[" << n << "]=" << p_nrm[n] << " " << p_numeric[n] << endl;
        EXPECT_NEAR(p_nrm[n],p_numeric[n],relative_error*p_numeric[n]);
    }
}
#endif

#endif // DO_THIS_TEST