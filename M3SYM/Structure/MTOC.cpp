
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