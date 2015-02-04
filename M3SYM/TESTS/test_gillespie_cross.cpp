
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

//#define DO_THIS_LONG_TEST

#ifdef DO_THIS_LONG_TEST 

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>

using namespace boost::accumulators;

#include "gtest/gtest.h"

#include "common.h"

#include "Species.h"
#include "Reaction.h"
#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"
#include "ChemSim.h"

vector<double> A1_A8_Network (int method)
{
    const long long int N_SAMPLE_POINTS=pow(10,5);
    const long long int NA1MAX = 50;
    const double tau_snapshot = 1.9; //seconds
    
    SpeciesBulk A1("A1", 35);
    SpeciesBulk A2("A2", 10);
    SpeciesBulk A3("A3",  8);
    SpeciesBulk A4("A4",  0);
    SpeciesBulk A5("A5",  3);
    SpeciesBulk A6("A6",  2);
    SpeciesBulk A7("A7",  1);
    SpeciesBulk A8("A8",  0);
    
    Reaction r1f = { {&A1,&A2, &A3,&A4}, 2, 2, 19.1}; // A1 + A2 -> A3 + A4
    Reaction r1b = { {&A3,&A4, &A1,&A2}, 2, 2, 23.3}; // A3 + A4 -> A1 + A2
    Reaction r2f = { {&A2, &A5}, 1, 1, 1.2};          // A2 -> A5
    Reaction r2b = { {&A5, &A2}, 1, 1, 3.2};          // A5 -> A2
    Reaction r3 = { {&A3, &A6, &A7}, 1, 2, 3.9};      // A3 -> A6 + A7
    Reaction r4 = { {&A6, &A1}, 1, 1, 10.9};          // A6 -> A1
    Reaction r5 = { {&A2, &A3}, 1, 1, 2.3};           // A2 -> A3
    Reaction r6 = { {&A5, &A6, &A7, &A8}, 2, 2, 3.9}; // A5 + A6 -> A7 + A8
    Reaction r7 = { {&A7, &A8, &A1}, 2, 1, 0.9};      // A7 + A8 -> A1
    Reaction r8 = { {&A2, &A7, &A8}, 1, 2, 8.9};      // A2 -> A7 + A8
    Reaction r9 = { {&A1, &A2}, 1, 1, 12.4};          // A1 -> A2
    Reaction r10 = { {&A4, &A1}, 1, 1, 16.4};         // A4 -> A1
    
    ChemSimImpl *chem_sim_impl = nullptr;
    
    switch (method) {
        case 0:
            chem_sim_impl = new ChemNRMImpl;
            break;
        case 1:
            chem_sim_impl = new ChemGillespieImpl;
            break;
        case 2:
            chem_sim_impl = new ChemSimpleGillespieImpl;
            break;
        default:
            assert(0 && "The method variable can only be 0, 1, or 2.");
    }
    ChemSim::setInstance(ChemSimInitKey(), chem_sim_impl);
    

    ChemSim::addReaction(ChemSimReactionKey(), &r3);
    ChemSim::addReaction(ChemSimReactionKey(), &r4);
    ChemSim::addReaction(ChemSimReactionKey(), &r5);
    ChemSim::addReaction(ChemSimReactionKey(), &r6);
    ChemSim::addReaction(ChemSimReactionKey(), &r7);
    ChemSim::addReaction(ChemSimReactionKey(), &r8);
    ChemSim::addReaction(ChemSimReactionKey(), &r9);
    ChemSim::addReaction(ChemSimReactionKey(), &r10);
    ChemSim::addReaction(ChemSimReactionKey(), &r1f);
    ChemSim::addReaction(ChemSimReactionKey(), &r1b);
    ChemSim::addReaction(ChemSimReactionKey(), &r2f);
    ChemSim::addReaction(ChemSimReactionKey(), &r2b);

    ChemSim::initialize(ChemSimInitKey());
    
    vector<long long int> x_hist(NA1MAX+1);
    
    for(long long int i=0;i<N_SAMPLE_POINTS;++i){
        A1.setN(35);
        A2.setN(10);
        A3.setN(8);
        A4.setN(0);
        A5.setN(3);
        A6.setN(2);
        A7.setN(1);
        A8.setN(0);
        
        long long int x_pentult=0;
        chem.initialize();
        long long int events=0;
        do {
            x_pentult=A1.getN();
            bool success = ChemSim::run(ChemSimRunKey(), 1);
            if(!success){
                cout << "chem.run(1) has failed, i= " << i << endl;
                ChemSim::printReactions();
                break;
            }
            ++events;
        } while (tau()<tau_snapshot);
        ++x_hist[x_pentult];
        if(i%10000==0)
            cout << "For loop, i=" << i << ", events=" << events << endl;
    }
    
    
    vector<double> p_nrm;
    
    for(int n=0; n<(NA1MAX+1); ++n){
        double p_est=double(x_hist[n])/N_SAMPLE_POINTS;
        p_nrm.push_back(p_est);
    }
    return p_nrm;
}


TEST(GillespieCrossTest, A1_A8_Network) {
    auto p_nrm = A1_A8_Network(0);
    auto p_gillespie = A1_A8_Network(1);
    auto p_simple_gillespie = A1_A8_Network(2);
    
    double relative_error=0.05; //i.e. allow a 5% relative error
    
    for(int n=4; n<16; ++n){
        EXPECT_NEAR(p_nrm[n],p_gillespie[n],relative_error*p_nrm[n]);
        EXPECT_NEAR(p_nrm[n],p_simple_gillespie[n],relative_error*p_simple_gillespie[n]);
    }
    
}



#endif

//SimpleGillespie (132 s)    NRM (220 s)       Caching Gillespie (122 s)
//
//P[0]=2e-05               P[0]=4e-05        P[0]=0
//P[1]=0.00021  	  	   P[1]=0.00022    	   P[1]=0.00021
//P[2]=0.0013   	  	   P[2]=0.00132    	   P[2]=0.00128
//P[3]=0.00475  	  	   P[3]=0.0047     	   P[3]=0.00488
//P[4]=0.01415  	  	   P[4]=0.01445    	   P[4]=0.01445
//P[5]=0.03132  	  	   P[5]=0.03223    	   P[5]=0.03078
//P[6]=0.05758  	  	   P[6]=0.05779    	   P[6]=0.05678
//P[7]=0.08984  	  	   P[7]=0.08991    	   P[7]=0.08892
//P[8]=0.11406  	  	   P[8]=0.11726    	   P[8]=0.11711
//P[9]=0.13629  	  	   P[9]=0.13618    	   P[9]=0.13673
//P[10]=0.1374  	  	   P[10]=0.13715   	   P[10]=0.13813
//P[11]=0.12517 	  	   P[11]=0.12447   	   P[11]=0.12352
//P[12]=0.10306 	  	   P[12]=0.10058   	   P[12]=0.10166
//P[13]=0.0743  	  	   P[13]=0.07422   	   P[13]=0.07366
//P[14]=0.0488  	  	   P[14]=0.04869   	   P[14]=0.04923
//P[15]=0.03005 	  	   P[15]=0.02973   	   P[15]=0.03062
//P[16]=0.01636 	  	   P[16]=0.01594   	   P[16]=0.01696
//P[17]=0.00886 	  	   P[17]=0.00867   	   P[17]=0.0083
//P[18]=0.00399 	  	   P[18]=0.00381   	   P[18]=0.00389
//P[19]=0.00153 	  	   P[19]=0.00166   	   P[19]=0.00178
//P[20]=0.0007  	  	   P[20]=0.00058   	   P[20]=0.00073
//P[21]=0.00017 	  	   P[21]=0.00026   	   P[21]=0.00023
//P[22]=8e-05   	  	   P[22]=0.00013   	   P[22]=0.00012
//P[23]=1e-05   	  	   P[23]=1e-05     	   P[23]=2e-05
//P[24]=0	      	  	   P[24]=0	   	   P[24]=1e-05
//P[25]=0	      	  	   P[25]=0         	   P[25]=0
//
