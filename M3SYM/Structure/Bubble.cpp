
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

#include "Bubble.h"

#include "SubSystem.h"
#include "Bead.h"

Bubble::Bubble(SubSystem* ps, vector<double> coordinates,
               short type, double radius, double kRepuls, double screenLength)

    : Trackable(true, false, true, false), _ps(ps), _type(type),
      _radius(radius), _kRepuls(kRepuls), _screenLength(screenLength),
      _ID(_bubbles.getID()), coordinate(coordinates) {
    
    //set up bead
    _bead = _ps->addTrackable<Bead>(coordinates, this, 0);
}

void Bubble::updatePosition() {
    
    coordinate = _bead->coordinate;
}

void Bubble::printSelf() {
    
    cout << endl;
    
    cout << "Bubble: ptr = " << this << endl;
    cout << "Bubble ID = " << _ID << endl;
    cout << "Bubble type = " << _type << endl;
    cout << "Bubble radius = " << _radius << endl;
    
    cout << endl;
    
    cout << "Bead information..." << endl;
    
    _bead->printSelf();
    
    cout << endl;
}

Database<Bubble*> Bubble::_bubbles;
