
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

#include "FilamentInitializer.h"

#include "Boundary.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"

using namespace mathfunc;

vector<tuple<short, vector<double>, vector<double>>>
RandomFilamentDist::createFilaments(Boundary* b, int numFilaments, short filamentType, int lenFilaments) {
    
    vector<tuple<short, vector<double>, vector<double>>> filamentData;
    
    //Create random distribution of filaments
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        double firstX = Rand::randDouble(0,1) * SysParams::Geometry().NX * SysParams::Geometry().compartmentSizeX;
        double firstY = Rand::randDouble(0,1) * SysParams::Geometry().NY * SysParams::Geometry().compartmentSizeY;
        double firstZ = Rand::randDouble(0,1) * SysParams::Geometry().NZ * SysParams::Geometry().compartmentSizeZ;
        
        double directionX = Rand::randDouble(-1,1);
        double directionY = Rand::randDouble(-1,1);
        double directionZ = Rand::randDouble(-1,1);
        
        //Create a random filament vector one cylinder long
        vector<double> firstPoint = {firstX, firstY, firstZ};
        
        vector<double> direction = {directionX, directionY, directionZ};
        normalize(direction);
        
        vector<double> secondPoint =
            nextPointProjection(firstPoint,(double)lenFilaments *
            SysParams::Geometry().cylinderSize[filamentType] - 0.01, direction);
        
        if(b->within(firstPoint) && b->within(secondPoint)) {
            filamentData.emplace_back(filamentType, firstPoint, secondPoint);
            filamentCounter++;
        }
    }
    return filamentData;
}

