
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
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
    
    vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> filaments;
    vector<tuple<string, short, vector<vector<floatingpoint>>>> dummy;
    vector<tuple<string, short, vector<floatingpoint>>> dummy2;
    vector<vector<floatingpoint>> dummy3;
    //Create random distribution of filaments
    int filamentCounter = 0;

    //if boundary shape is cylinder, create filament in the center of system and perpendicular to Z axis
    if(b->getShape() == BoundaryShape::Cylinder) {

        while (filamentCounter < numFilaments) {

            //Create a random filament vector one cylinder long
            vector<floatingpoint> firstPoint = GController::getRandomCenterCoordinates();

            floatingpoint directionX = Rand::randfloatingpoint(-1,1);
            floatingpoint directionY = Rand::randfloatingpoint(-1,1);

            floatingpoint directionZ = 0;
            vector<floatingpoint> direction = normalizeVector({directionX, directionY, directionZ});

            vector<floatingpoint> secondPoint =
                nextPointProjection(firstPoint,(floatingpoint)lenFilaments *
                SysParams::Geometry().cylinderSize[filamentType] - 0.01, direction);

            //check if these points are outside bubbles
            bool inBubble = false;
            for(auto bb : Bubble::getBubbles()) {
                auto radius = bb->getRadius();

                if((twoPointDistancesquared(bb->getBead()->vcoordinate(), firstPoint) < (radius * radius)) ||
                   (twoPointDistancesquared(bb->getBead()->vcoordinate(), secondPoint) < (radius * radius)))
                    inBubble = true;
            }

            //check if within cutoff of boundary
            bool outsideCutoff = false;
            if(b->distance(firstPoint) < SysParams::Boundaries().BoundaryCutoff / 4.0 ||
               b->distance(secondPoint) < SysParams::Boundaries().BoundaryCutoff / 4.0) {
                outsideCutoff = true;
            }

            if(b->within(firstPoint) && b->within(secondPoint) && !inBubble && !outsideCutoff) {
                filaments.emplace_back(filamentType, firstPoint, secondPoint);
                filamentCounter++;
            }
        }
        return make_tuple(filaments, dummy, dummy2, dummy3);
    }

    else{
        while (filamentCounter < numFilaments) {

            //Create a random filament vector one cylinder long
            vector<floatingpoint> firstPoint = GController::getRandomCoordinates();

            floatingpoint directionX = Rand::randfloatingpoint(-1,1);
            floatingpoint directionY = Rand::randfloatingpoint(-1,1);
            floatingpoint directionZ = Rand::randfloatingpoint(-1,1);
            vector<floatingpoint> direction = normalizeVector({directionX, directionY, directionZ});

            vector<floatingpoint> secondPoint =
            nextPointProjection(firstPoint,(floatingpoint)lenFilaments *
                                SysParams::Geometry().cylinderSize[filamentType] - 0.01,
                                direction);

            //check if these points are outside bubbles
            bool inBubble = false;
            for(auto bb : Bubble::getBubbles()) {
                auto radius = bb->getRadius();
                
                if((twoPointDistancesquared(bb->getBead()->vcoordinate(), firstPoint) < (radius * radius)) ||
                   (twoPointDistancesquared(bb->getBead()->vcoordinate(), secondPoint) < (radius * radius)))
                    inBubble = true;
            }

            //check if within cutoff of boundary
            bool outsideCutoff = false;
            if(b->distance(firstPoint) < SysParams::Boundaries().BoundaryCutoff / 4.0 ||
               b->distance(secondPoint) < SysParams::Boundaries().BoundaryCutoff / 4.0) {
                outsideCutoff = true;
            }

            if(b->within(firstPoint) && b->within(secondPoint) && !inBubble && !outsideCutoff) {
                filaments.emplace_back(filamentType, firstPoint, secondPoint);
                filamentCounter++;
            }
        }
        return make_tuple(filaments, dummy, dummy2, dummy3);

    }



}
FilamentData ConnectedFilamentDist::createFilaments(Boundary* b, int numFilaments,
                                                    int filamentType,
                                                    int lenFilaments) {

    ///SET THIS SPACING PARAMETER
    floatingpoint maxSpacing = 50;
    
    vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> filaments;
    vector<tuple<string, short, vector<vector<floatingpoint>>>> dummy;
    vector<tuple<string, short, vector<floatingpoint>>> dummy2;
    vector<vector<floatingpoint>> dummy3;
    
    ///First filament as normal
    //Create a random filament vector one cylinder long
    vector<floatingpoint> firstPoint = {500,1000,1000};
    
    vector<floatingpoint> direction = normalizeVector({1, 0, 0});
    
    vector<floatingpoint> secondPoint =
    nextPointProjection(firstPoint,(floatingpoint)lenFilaments *
                        SysParams::Geometry().cylinderSize[filamentType] - 0.01, direction);
    
    floatingpoint len = twoPointDistance(firstPoint, secondPoint);
    filaments.emplace_back(filamentType, firstPoint, secondPoint);
    
    vector<floatingpoint> prevFirstPoint = firstPoint;
    vector<floatingpoint> prevSecondPoint = secondPoint;
    
    floatingpoint safeDist = SysParams::Boundaries().BoundaryCutoff;
    
    ///now create properly distanced network
    int filamentCounter = 1;
    while (filamentCounter < numFilaments) {
    
        ///pick a random distance from a random point on the chain
        direction = twoPointDirection(firstPoint, secondPoint);
        len = twoPointDistance(firstPoint, secondPoint);
        floatingpoint randomSeg = Rand::randfloatingpoint(0, len);
        
        vector<floatingpoint> randomPoint = nextPointProjection(firstPoint, randomSeg, direction);
        
        //now pick another random point which is within a certain distance away
        floatingpoint directionX = Rand::randfloatingpoint(-1,1);
        floatingpoint directionY = Rand::randfloatingpoint(-1,1);
        floatingpoint directionZ = Rand::randfloatingpoint(-1,1);
        vector<floatingpoint> randDirection = normalizeVector({directionX, directionY, directionZ});
        
        floatingpoint randomDist = Rand::randfloatingpoint(0, maxSpacing);
        vector<floatingpoint> nextRandomPoint = nextPointProjection(randomPoint, randomDist, randDirection);
        
        //now pick another random direction for the next filament creation
        directionX = Rand::randfloatingpoint(-1,1);
        directionY = Rand::randfloatingpoint(-1,1);
        directionZ = Rand::randfloatingpoint(-1,1);
        randDirection = normalizeVector({directionX, directionY, directionZ});
    
        //random length spacing for new filament
        floatingpoint randomLengthSpacing = Rand::randfloatingpoint(0, (floatingpoint)lenFilaments *
                                     SysParams::Geometry().cylinderSize[filamentType] - 0.01);
        
        firstPoint = nextPointProjection(nextRandomPoint, randomLengthSpacing, randDirection);
        
        //switch rand direction and create second point
        randDirection[0] = - randDirection[0];
        randDirection[1] = - randDirection[1];
        randDirection[2] = - randDirection[2];
        
        secondPoint = nextPointProjection(firstPoint, (floatingpoint)lenFilaments *
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
    
    return make_tuple(filaments,dummy,dummy2, dummy3);
}

FilamentData MTOCFilamentDist::createFilaments(Boundary* b, int numFilaments,
                                                            int filamentType,
                                                            int lenFilaments) {
    
    vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> filaments;
    vector<tuple<string, short, vector<vector<floatingpoint>>>> dummy;
    vector<vector<floatingpoint>> dummy3;
    vector<tuple<string, short, vector<floatingpoint>>> dummy2;
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        floatingpoint l = Rand::randfloatingpoint(0,2 * M_PI);
        floatingpoint h = Rand::randfloatingpoint(-M_PI/2, M_PI/2);
        
        vector<floatingpoint> point1;
        point1.push_back(_coordMTOC[0] + _radius * cos(l) * cos(h));
        point1.push_back(_coordMTOC[1] + _radius * sin(h));
        point1.push_back(_coordMTOC[2] + _radius * sin(l) * cos(h));
        
        // get projection outward from the MTOC
        auto dir = normalizeVector(twoPointDirection(_coordMTOC, point1));
        auto point2 = nextPointProjection(point1,
            SysParams::Geometry().cylinderSize[filamentType]*lenFilaments - 0.01, dir);
        
        filaments.emplace_back(filamentType, point1, point2);
        
        filamentCounter++;
    }
    
    return make_tuple(filaments,dummy,dummy2, dummy3);
}
