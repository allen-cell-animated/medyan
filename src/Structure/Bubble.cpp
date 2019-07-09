
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

#include "Bubble.h"

#include "SubSystem.h"
#include "Bead.h"

#include "SysParams.h"

Bubble::Bubble(SubSystem* ps, vector<floatingpoint> coordinates, short type)

    : Trackable(true, false, true, false), _ps(ps), _type(type),
      coordinate(coordinates) {
    
    //set up mechanical constants
    _kRepuls = SysParams::Mechanics().BubbleK[_type];
    _radius  = SysParams::Mechanics().BubbleRadius[_type];
    _screenLength = SysParams::Mechanics().BubbleScreenLength[_type];
    
    if(SysParams::Mechanics().MTOCBendingK.size() == 0){
        _MTOCBendingK = 0.0;
    }
    else{
      _MTOCBendingK = SysParams::Mechanics().MTOCBendingK[_type];
    }
          
    //set up bead
    _bead = _ps->addTrackable<Bead>(coordinates, this, 0);
}

void Bubble::updatePosition() {
    
    coordinate = _bead->vcoordinate();
}

void Bubble::printSelf()const {
    
    cout << endl;
    
    cout << "Bubble: ptr = " << this << endl;
    cout << "Bubble ID = " << getId() << endl;
    cout << "Bubble type = " << _type << endl;
    cout << "Bubble radius = " << _radius << endl;
    
    cout << endl;
    
    cout << "Bead information..." << endl;
    
    _bead->printSelf();
    
    cout << endl;
}
