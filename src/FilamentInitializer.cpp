
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
#include "MembraneRegion.h"

#include "MathFunctions.h"
#include "GController.h"
#include "SysParams.h"
#include "Rand.h"

using namespace mathfunc;

FilamentData RandomFilamentDist::createFilaments(const MembraneRegion& mr,
                                                 int numFilaments,
                                                 int filamentType,
                                                 int lenFilaments) {
    
    vector<tuple<short, vector<double>, vector<double>>> filaments;
    vector<tuple<string, short, vector<vector<double>>>> dummy;
    vector<tuple<string, short, vector<double>>> dummy2;
    vector<vector<double>> dummy3;
    //Create random distribution of filaments
    int filamentCounter = 0;
    
    Boundary *b = mr.getBoundary();
    
    //Qin, if boundary shape is cylinder, create filament in the center of system and vertical to Z axis
    if(b && b->getShape() == BoundaryShape::Cylinder) {
        
        while (filamentCounter < numFilaments) {

            //Create a random filament vector one cylinder long
            vector<double> firstPoint = GController::getRandomCenterCoordinates();

            double directionX = Rand::randDouble(-1,1);
            double directionY = Rand::randDouble(-1,1);

            double directionZ = 0;
            vector<double> direction = normalizeVector({directionX, directionY, directionZ});

            vector<double> secondPoint =
                nextPointProjection(firstPoint,(double)lenFilaments *
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

    //Qin
    else{
        while (filamentCounter < numFilaments) {

            //Create a random filament vector one cylinder long
            vector<double> firstPoint = GController::getRandomCoordinates();

            double directionX = Rand::randDouble(-1,1);
            double directionY = Rand::randDouble(-1,1);
            double directionZ = Rand::randDouble(-1,1);
            vector<double> direction = normalizeVector({directionX, directionY, directionZ});

            vector<double> secondPoint =
            nextPointProjection(firstPoint,(double)lenFilaments *
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
            bool outsideCutoff = mr.getBoundary() && (
                mr.getBoundary()->distance(firstPoint) < SysParams::Boundaries().BoundaryCutoff / 4.0 ||
                mr.getBoundary()->distance(secondPoint) < SysParams::Boundaries().BoundaryCutoff / 4.0
            );
            
            if(mr.contains(vector2Vec<3, double>(firstPoint)) && mr.contains(vector2Vec<3, double>(secondPoint)) && !inBubble && !outsideCutoff) {
                filaments.emplace_back(filamentType, firstPoint, secondPoint);
                filamentCounter++;
            }
        }
        return make_tuple(filaments, dummy, dummy2, dummy3);

    }



}

FilamentData ConnectedFilamentDist::createFilaments(const MembraneRegion& mr,
                                                    int numFilaments,
                                                    int filamentType,
                                                    int lenFilaments) {

    ///SET THIS SPACING PARAMETER
    double maxSpacing = 50;
    
    vector<tuple<short, vector<double>, vector<double>>> filaments;
    vector<tuple<string, short, vector<vector<double>>>> dummy;
    vector<tuple<string, short, vector<double>>> dummy2;
    vector<vector<double>> dummy3;
    
    ///First filament as normal
    //Create a random filament vector one cylinder long
    vector<double> firstPoint = {500,1000,1000};
    
    vector<double> direction = normalizeVector({1, 0, 0});
    
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
        vector<double> randDirection = normalizeVector({directionX, directionY, directionZ});
        
        double randomDist = Rand::randDouble(0, maxSpacing);
        vector<double> nextRandomPoint = nextPointProjection(randomPoint, randomDist, randDirection);
        
        //now pick another random direction for the next filament creation
        directionX = Rand::randDouble(-1,1);
        directionY = Rand::randDouble(-1,1);
        directionZ = Rand::randDouble(-1,1);
        randDirection = normalizeVector({directionX, directionY, directionZ});
    
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
        if(
            mr.contains(vector2Vec<3, double>(firstPoint)) && mr.contains(vector2Vec<3, double>(secondPoint)) && (
                !mr.getBoundary() || (
                    mr.getBoundary()->distance(firstPoint) > safeDist &&
                    mr.getBoundary()->distance(secondPoint) > safeDist
                )
            )
        ) {
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

FilamentData MTOCFilamentDist::createFilaments(const MembraneRegion& mr,
                                               int numFilaments,
                                               int filamentType,
                                               int lenFilaments) {
    
    vector<tuple<short, vector<double>, vector<double>>> filaments;
    vector<tuple<string, short, vector<vector<double>>>> dummy;
    vector<vector<double>> dummy3;
    vector<tuple<string, short, vector<double>>> dummy2;
    int filamentCounter = 0;
    while (filamentCounter < numFilaments) {
        
        double l = Rand::randDouble(0,2 * M_PI);
        double h = Rand::randDouble(-M_PI/2, M_PI/2);
        
        vector<double> point1;
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
