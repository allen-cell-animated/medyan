
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CGPolakRibiereMethod_h
#define MEDYAN_CGPolakRibiereMethod_h

#include <cmath>
#include <numeric>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"

/// The Polak-Ribiere method for conjugate gradient minimization
class PolakRibiere : public CGMethod
{
public:
    
    virtual MinimizationResult minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                          floatingpoint MAXDIST, floatingpoint LAMBDAMAX,
                          floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY, string _LINESEARCHALGORITHM,
                          bool steplimit);
protected:
	chrono::high_resolution_clock::time_point tbegin, tend;
	floatingpoint prevlambda = 0;
	floatingpoint prevbeta = 0;
	unsigned int skipcounter = 0;

	void calculateEvsalpha(ForceFieldManager &FFM, floatingpoint lambda, floatingpoint
	LAMBDAMAX, floatingpoint FDotFA);
#ifdef CUDAACCL
    cudaStream_t stream_shiftsafe = NULL, stream_dotcopy = NULL;
    cudaStream_t stream1 = NULL, stream2 = NULL, stream3 = NULL;
    cudaEvent_t  event_safe = NULL, event_dot = NULL;
    cudaStream_t Ms1 = NULL, Ms2 = NULL, Ms3 = NULL, Ms4 = NULL, *Msp1 = NULL, *Msp2 = NULL, *Msps = NULL;
    cudaEvent_t Me1 = NULL, Me2 = NULL, *Mep1 = NULL, *Mep2 = NULL, *Meps = NULL;
#endif
};
#endif

