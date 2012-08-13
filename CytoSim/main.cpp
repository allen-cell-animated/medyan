//
//  main.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

/*! \mainpage The Main Page for the CytoSim software package 
 
 \section intro_sec Introduction
 
 Cell motility plays a key role in human biology and disease, contributing ubiquitously to such important processes as embryonic development, wound repair and cancer metastasis. Papoian laboratory is interested in gaining deeper understanding of the physical chemistry behind these complex, far-from-equilibrium mechano-chemical processes. His approach is based on combining stochastic reaction-diffusion treatment of cellular biochemical processes with polymer physics of cytoskeletal filament network growth, while explicitly coupling chemistry and mechanics. 
 
 Papoian laboratory has developed **CytoSim**, a software package to simulate growth dynamics of actin based filamentous networks *in vitro* and *in vivo*. Recent papers where **CytoSim** or its predecessor, **StochTools**, were used can be found on the publication section of [the Papoian group's main web page: ](http://papoian.chem.umd.edu/)
 
 
 \section install_sec Installation
 
 \subsection step1 Step 1: Prerequisites
 
 The following software packages need to be installed first:
 
 - Boost 1.49
 - GSL ...
 
 \subsection step2 Step 2: Installation of CytoSim itself
 
 Untar the **CytoSim** source code into some directory, enter into the "CytoSim/CytoSim" and execute "make" from the command line.
 
 
 */

#include <iostream>

#include "Species.h"
#include "Reaction.h"
#include "ChemGillespieImpl.h"
//#include "ChemNRMImpl.h"
#include "ChemSim.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>

using namespace boost::accumulators;


using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    const long long int N_SAMPLE_POINTS=pow(10,6);
    const long long int Nstart = 10;
    const double tau_snapshot = 0.48; //seconds
    //long long int print_freq = pow(10,7);
    SpeciesBulk A1("A1",  Nstart);
    SpeciesBulk A2("A2", 0);
    // A1 <-> A2 with the same forward and backward rates; [A]~[B] at steady state
    Reaction r1 = { {&A1,&A2}, 1, 1, 2.5 };
    Reaction r2 = { {&A2,&A1}, 1, 1, 2.5 };
    
    ChemGillespieImpl chem_impl;
    ChemSim chem(&chem_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    //chem.printReactions();
    
    r1.printDependents();
    cout << endl << endl;
    r2.printDependents();
//    return 0;
    
    vector<long long int> n_hist(Nstart+1);
    
    //int counter = 0;
    accumulator_set<double, stats<tag::mean>> accTau;
    long long int N_penultimate;
    for(long long int i=0;i<N_SAMPLE_POINTS;++i){
        A1.setN(Nstart);
        A2.setN(0);
        chem.initialize();
        //        long long int events=0;
        if(i%100000==0)
            cout << i << endl;
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
        cout << "P[" << n << "]=" << p_est << " " << p_analyt << endl;
        // EXPECT_NEAR(p_est,p_analyt,relative_error*p_analyt);
    }
    // Note that the error for P[0] is hard to converge, most likely because of P[0] being very small
    
    cout << "Main exited..." << endl;
    return 0;
}

