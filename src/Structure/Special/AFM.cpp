
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "AFM.h"

#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

void AFM::setBubble(Bubble* b) {
    
    _bubble = b;
    _bubble->setAsAFM();
    
    addChild(unique_ptr<Component>(b));
}

void AFM::printSelf() const {
    
    cout << endl;
    
    cout << "AFM: ptr = " << this << endl;
    cout << "AFM ID = " << getId() << endl;
    
    cout << endl;
    cout << "Bubble information..." << endl;
    
    _bubble->printSelf();
    
    cout << endl;
}
