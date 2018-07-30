
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "FilamentInitializer.h"

#include "Boundary.h"
#include "Bead.h"

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

FilamentData ConnectedFilamentDist::createFilaments(Boundary* b, int numFilaments,
                                                    int filamentType,
                                                    int lenFilaments) {

    ///SET THIS SPACING PARAMETER
    double maxSpacing = 50;
    
    FilamentData filaments;
    
    ///First filament as normal
    //Create a random filament vector one cylinder long
    vector<double> firstPoint = {500,1000,1000};
    
    vector<double> direction = normalizedVector({1, 0, 0});
    
    vector<double> secondPoint =
    nextPointProjection(firstPoint,(double)lenFilaments *
                        SysParams::Geometry().cylinderSize[filamentType] - 0.01, direction);
    
    double len = twoPointDistance(firstPoint, secondPoint);
    filaments.emplace_back(filamentType, firstPoint, secondPoint);
    
    vector<double> prevFirstPoint = firstPoint;
    vector<double> prevSecondPoint = secondPoint;
    
    double safeDist = SysParams::Boundaries().BoundaryCutoff;
    
    ///now create properly distanced network
    int filamentCounter = 1;
    while (filamentCounter < numFilaments) {
    
        ///pick a random distance from a random point on the chain
        direction = twoPointDirection(firstPoint, secondPoint);
        len = twoPointDistance(firstPoint, secondPoint);
        double randomSeg = Rand::randDouble(0, len);
        
        vector<double> randomPoint = nextPointProjection(firstPoint, randomSeg, direction);
        
        //now pick another random point which is within a certain distance away
        double directionX = Rand::randDouble(-1,1);
        double directionY = Rand::randDouble(-1,1);
        double directionZ = Rand::randDouble(-1,1);
        vector<double> randDirection = normalizedVector({directionX, directionY, directionZ});
        
        double randomDist = Rand::randDouble(0, maxSpacing);
        vector<double> nextRandomPoint = nextPointProjection(randomPoint, randomDist, randDirection);
        
        //now pick another random direction for the next filament creation
        directionX = Rand::randDouble(-1,1);
        directionY = Rand::randDouble(-1,1);
        directionZ = Rand::randDouble(-1,1);
        randDirection = normalizedVector({directionX, directionY, directionZ});
    
        //random length spacing for new filament
        double randomLengthSpacing = Rand::randDouble(0, (double)lenFilaments *
                                     SysParams::Geometry().cylinderSize[filamentType] - 0.01);
        
        firstPoint = nextPointProjection(nextRandomPoint, randomLengthSpacing, randDirection);
        
        //switch rand direction and create second point
        randDirection[0] = - randDirection[0];
        randDirection[1] = - randDirection[1];
        randDirection[2] = - randDirection[2];
        
        secondPoint = nextPointProjection(firstPoint, (double)lenFilaments *
                      SysParams::Geometry().cylinderSize[filamentType] - 0.01, randDirection);
        
        //choose if within boundary
        if(b->within(firstPoint) && b->within(secondPoint) &&
           b->distance(firstPoint) > safeDist &&
           b->distance(secondPoint) > safeDist) {
            filaments.emplace_back(filamentType, firstPoint, secondPoint);
            
            prevFirstPoint = firstPoint;
            prevSecondPoint = secondPoint;
            
            filamentCounter++;
        }
        else { //reset
            firstPoint = prevFirstPoint;
            secondPoint = prevSecondPoint;
        }
        
    }
    
    return filaments;
}
