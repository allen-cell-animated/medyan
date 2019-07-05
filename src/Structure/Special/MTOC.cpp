
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MTOC.h"

#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

void MTOC::setBubble(Bubble* b) {
    
    _bubble = b;
    _bubble->setAsMTOC();
    
    addChild(unique_ptr<Component>(b));
}

void MTOC::printSelf() {
    
    cout << endl;
    
    cout << "MTOC: ptr = " << this << endl;
    cout << "MTOC ID = " << _ID << endl;
    
    cout << endl;
    cout << "Bubble information..." << endl;
    
    _bubble->printSelf();
    
    cout << endl;
}


Database<MTOC*> MTOC::_mtocs;
