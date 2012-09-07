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
//#define TRACK_UPPER_COPY_N

#include <cstdint>

typedef unsigned short species_copy_t;
const species_copy_t max_ulim = 1024;
extern double global_time;

inline double tau() {return global_time;}

#endif
