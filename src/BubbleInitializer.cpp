
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

#include "BubbleInitializer.h"

#include "Boundary.h"

#include "GController.h"
#include "SysParams.h"
#include "Rand.h"

BubbleData RandomBubbleDist::createBubbles(Boundary* b, int numBubbles,
                                                        int bubbleType) {
    
    BubbleData bubbles;
    
    //Create random distribution of Bubbles
    int bubbleCounter = 0;
    while (bubbleCounter < numBubbles) {
        
        //Create a random Bubble coord
        vector<floatingpoint> coord = GController::getRandomCoordinates();
        
        if(b->within(coord) &&
           b->distance(coord) > SysParams::Mechanics().BubbleRadius[bubbleType]) {
            
            bubbles.emplace_back(bubbleType, coord);
            bubbleCounter++;
        }
    }
    return bubbles;
}
