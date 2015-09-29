
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

#include "common.h"

double global_time;

vector<short> B_BINDING_INDEX = vector<short>(MAX_FILAMENTS);
vector<short> L_BINDING_INDEX = vector<short>(MAX_FILAMENTS);
vector<short> M_BINDING_INDEX = vector<short>(MAX_FILAMENTS);

vector<vector<short>> BINDING_INDEX = vector<vector<short>>(MAX_FILAMENTS);