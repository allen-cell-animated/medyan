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
//#include "ChemSimpleGillespieImpl.h"
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
    const long long int N_SAMPLE_POINTS=pow(10,5);
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
    Reaction xa = { {&X,&A}, 1, 1, kxa };
    Reaction ax = { {&A,&X}, 1, 1, kax };
    Reaction r1 = { {&A,&B}, 1, 1, kab };
    Reaction r2 = { {&B,&C}, 1, 1, kbc };
    Reaction r3 = { {&C,&A}, 1, 1, kca };
    
    ChemGillespieImpl chem_nrm_impl;
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
            //            cout << "tau=" << tau() << ", X=" << X.getN() << ", A=" << A.getN() << ", B=" << B.getN() << ", C=" << C.getN() << endl;
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
    //    double relative_error=0.05; //i.e. allow a 5% relative error
    
    for(int n=0; n<(Nstart+4); ++n){
        cout << "P[" << n << "]=" << p_nrm[n] << " " << p_numeric[n] << endl;
    }
    
    cout << "Main exited..." << endl;
    return 0;
}

