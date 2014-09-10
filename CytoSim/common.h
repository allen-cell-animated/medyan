//
//  common.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_common_h
#define CytoSim_Experimenting_common_h

#define TRACK_DEPENDENTS
#define TRACK_ZERO_COPY_N
#define TRACK_UPPER_COPY_N
#define RSPECIES_SIGNALING
#define REACTION_SIGNALING
//#define BOOST_MEM_POOL
//#define BOOL_POOL_NSIZE 65536

#include <cstdint>
#include <vector>
#include <iostream>
#include <math.h>
#include "MathFunctions.h"
#include "SystemParameters.h"


///species constants
typedef unsigned short species_copy_t;
const species_copy_t max_ulim = 10000;
extern double global_time;

inline double tau() {return global_time;}


#endif
