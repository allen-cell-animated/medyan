//
//  Mcommon.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_Mcommon_h
#define CytoMech_Mcommon_h

#include <iostream>
#include <math.h>
#include "MathFunctions.h"

class Bead;
class Filament;
class FilamentDB;
class BeadDB;
class SubSystem;
class Network;
class Cylinder;
class CylinderDB;

void FletcherRievesMethod(System* ps);
void PolakRibiereMethod(System* ps);

/// Some constants for potentials. Maybe not the best place for them, can be moved to other places;
const double L = 20.0;
const double teta = 180.0;
const double kS = 10.0;
const double kB = 10.0;
const double kTw = 10.0;


#endif
