
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
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
    
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                          double MAXDIST, double LAMBDAMAX, bool steplimit);

    
protected:
#ifdef CUDAACCL
    cudaStream_t stream_shiftsafe = NULL, stream_dotcopy = NULL;
    cudaStream_t stream1 = NULL, stream2 = NULL, stream3 = NULL;
    cudaEvent_t  event_safe = NULL, event_dot = NULL;
    cudaStream_t Ms1 = NULL, Ms2 = NULL, Ms3 = NULL, Ms4 = NULL, *Msp1 = NULL, *Msp2 = NULL, *Msps = NULL;
    cudaEvent_t Me1 = NULL, Me2 = NULL, *Mep1 = NULL, *Mep2 = NULL, *Meps = NULL;
#endif
};
#endif

