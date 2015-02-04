
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
#include "SystemParameters.h"

using namespace mathfunc;

vector<vector<vector<double>>>
RandomFilamentDist::createFilaments(Boundary* b, int numFilaments, int lenFilaments) {
    
    vector<vector<vector<double>>> filamentData;
    
    //Create random distribution of filaments
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        double firstX = randomDouble(0,1) *
            SystemParameters::Geometry().compartmentSizeX *
            SystemParameters::Geometry().NX;
        
        double firstY = randomDouble(0,1) *
            SystemParameters::Geometry().compartmentSizeY *
        SystemParameters::Geometry().NY;
        
        double firstZ = randomDouble(0,1) *
            SystemParameters::Geometry().compartmentSizeZ *
            SystemParameters::Geometry().NZ;
        
        double directionX = randomDouble(-1,1);
        double directionY = randomDouble(-1,1);
        double directionZ = randomDouble(-1,1);
        
        //Create a random filament vector one cylinder long
        vector<double> firstPoint = {firstX, firstY, firstZ};
        
        auto normFactor = sqrt(directionX * directionX +
                               directionY * directionY +
                               directionZ * directionZ);
        
        vector<double> direction = {directionX/normFactor,
                                    directionY/normFactor,
                                    directionZ/normFactor};
        
        vector<double> secondPoint =
            nextPointProjection(firstPoint,(double)lenFilaments *
            SystemParameters::Geometry().cylinderSize - 0.01, direction);
        
        if(b->within(firstPoint) && b->within(secondPoint)) {
            filamentData.push_back({firstPoint, secondPoint});
            filamentCounter++;
        }
    }
    return filamentData;
}

