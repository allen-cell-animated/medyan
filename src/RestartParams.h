
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

#ifndef MEDYAN_RestartParams_h
#define MEDYAN_RestartParams_h
#include <vector>
#include "utility.h"
#include "Cylinder.h"


using namespace std;
/*Structures to store data parsed from Restartfile*/
struct restartBeadData{
	vector<unsigned int> bsidvec; //stableid
	vector<int> filidvec;
	vector<int> filpos;
	vector<float> coordvec;
	vector<float> forceAuxvec;
};

/*Stores CylData for each cylinder */
struct restartCylData{
	unsigned int cylsid; //stableid
	unsigned int filid;
	unsigned int filtype;
	unsigned int filpos;
	unsigned int beadsidpairvec[2];
	unsigned int endstatusvec[2];
	unsigned int endmonomerpos[2];
	unsigned int totalmonomers;
	floatingpoint eqlen;
	Cylinder* cylinderpointer;
};

/* Stores FilData for each Filament */
struct restartFilData{
	unsigned int filid;
	unsigned int filType;
	vector<unsigned int> cylsidvec;
	Filament* filamentpointer;
};

#endif
