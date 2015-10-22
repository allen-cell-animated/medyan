
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
#include "GController.h"
#include "SysParams.h"
#include "Rand.h"

using namespace mathfunc;

FilamentData RandomFilamentDist::createFilaments(Boundary* b, int numFilaments,
                                                              int filamentType,
                                                              int lenFilaments) {
    
    FilamentData filaments;
    
    //Create random distribution of filaments
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        //Create a random filament vector one cylinder long
        vector<double> firstPoint = GController::getRandomCoordinates();
        
        double directionX = Rand::randDouble(-1,1);
        double directionY = Rand::randDouble(-1,1);
        double directionZ = Rand::randDouble(-1,1);
        vector<double> direction = normalizedVector({directionX, directionY, directionZ});
        
        vector<double> secondPoint =
            nextPointProjection(firstPoint,(double)lenFilaments *
            SysParams::Geometry().cylinderSize[filamentType] - 0.01, direction);
        
        if(b->within(firstPoint) && b->within(secondPoint)) {
            filaments.emplace_back(filamentType, firstPoint, secondPoint);
            filamentCounter++;
        }
    }
    return filaments;
}

