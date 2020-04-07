
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MotorGhostInteractions.h"

MotorGhost* MotorGhostInteractions::_motorCulprit;

vector<floatingpoint> MotorGhostInteractions::individualenergies;
vector<floatingpoint> MotorGhostInteractions::tpdistvec;
vector<floatingpoint> MotorGhostInteractions::eqlvec;
vector<floatingpoint> MotorGhostInteractions::kstrvec;