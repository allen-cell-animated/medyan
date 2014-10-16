//
//  common.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_common_h
#define CytoSim_Experimenting_common_h

///Will eventually be command line macros
//#define MECHANICS
#define CHEMISTRY

///If compiling for testing
#define TESTING

///Other chemistry macros
#define TRACK_DEPENDENTS
#define TRACK_ZERO_COPY_N
#define TRACK_UPPER_COPY_N
#define RSPECIES_SIGNALING
#define REACTION_SIGNALING
//#define BOOST_MEM_POOL
//#define BOOL_POOL_NSIZE 65536

///Some constants for optimization of monomer containers
#define NUMSPECIESFILAMENT          1
#define NUMSPECIESPLUSEND           1
#define NUMSPECIESMINUSEND          1
#define NUMSPECIESBOUND             1
#define NUMSPECIESLINKER            1
#define NUMSPECIESMOTOR             1

///species constants
typedef unsigned short species_copy_t;
const species_copy_t max_ulim = 10000;
extern double global_time;

inline double tau() {return global_time;}


#endif
