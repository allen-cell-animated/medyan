
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryElement.h"

Database<BoundaryElement*> BoundaryElement::_boundaryElements;


void BoundaryElement::printSelf() {
    
    cout << endl;
    
    cout << "Boundary element: ptr = " << this << endl;
    cout << "Coordinates = " << _coords[0] << ", " << _coords[1] << ", " << _coords[2] << endl;
    
    cout << endl;
}