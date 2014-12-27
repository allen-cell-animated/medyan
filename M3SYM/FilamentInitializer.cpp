
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include <random>

#include "FilamentInitializer.h"

#include "Boundary.h"

#include "MathFunctions.h"
#include "SystemParameters.h"

using namespace mathfunc;

vector<vector<vector<double>>>
RandomFilamentDist::createFilaments(Boundary* b, int numFilaments, int lenFilaments) {
    
    vector<vector<vector<double>>> filamentData;
    //Create random distribution of filaments
    default_random_engine generator;
    uniform_real_distribution<double> dU(0,1);
    uniform_real_distribution<double> dUNeg(-1,1);
    
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        double firstX = dU(generator) *
            SystemParameters::Geometry().compartmentSizeX *
            SystemParameters::Geometry().NX;
        double firstY = dU(generator) *
            SystemParameters::Geometry().compartmentSizeY *
        SystemParameters::Geometry().NY;
        double firstZ = dU(generator) *
            SystemParameters::Geometry().compartmentSizeZ *
            SystemParameters::Geometry().NZ;
        
        double directionX = dUNeg(generator);
        double directionY = dUNeg(generator);
        double directionZ = dUNeg(generator);
        
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

