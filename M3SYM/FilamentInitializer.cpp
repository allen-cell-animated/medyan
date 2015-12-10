
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
#include "Bubble.h"
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
        
        //check if these points are outside bubbles
        bool inBubble = false;
        for(auto bb : Bubble::getBubbles()) {
            
            if((twoPointDistance(bb->getBead()->coordinate, firstPoint) < bb->getRadius()) ||
               (twoPointDistance(bb->getBead()->coordinate, secondPoint) < bb->getRadius()))
                inBubble = true;
        }
        double safeDist = SysParams::Boundaries().BoundaryCutoff;
        
        if(b->within(firstPoint) && b->within(secondPoint) &&
           b->distance(firstPoint) > safeDist &&
           b->distance(secondPoint) > safeDist && !inBubble) {
            filaments.emplace_back(filamentType, firstPoint, secondPoint);
            filamentCounter++;
        }
    }
    return filaments;
}

FilamentData MTOCFilamentDist::createFilaments(Boundary* b, int numFilaments,
                                                            int filamentType,
                                                            int lenFilaments) {
    FilamentData filaments;
    
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        double l = Rand::randDouble(0,2 * M_PI);
        double h = Rand::randDouble(-M_PI/2, M_PI/2);
        
        vector<double> point1;
        point1.push_back(_coordMTOC[0] + _radius * cos(l) * cos(h));
        point1.push_back(_coordMTOC[1] + _radius * sin(h));
        point1.push_back(_coordMTOC[2] + _radius * sin(l) * cos(h));
        
        // get projection outward from the MTOC
        auto dir = normalizedVector(twoPointDirection(_coordMTOC, point1));
        auto point2 = nextPointProjection(point1,
            SysParams::Geometry().cylinderSize[filamentType]*lenFilaments - 0.01, dir);
        
        filaments.emplace_back(filamentType, point1, point2);
        
        filamentCounter++;
    }
    
    return filaments;
}
