
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

#include "Bubble.h"

#include "SubSystem.h"
#include "Bead.h"

#include "SysParams.h"
#include "CUDAcommon.h"

Bubble::Bubble(SubSystem* ps, vector<double> coordinates, short type)

    : Trackable(true, false, true, false), _ps(ps), _type(type),
      _ID(_bubbles.getID()), coordinate(coordinates) {
    
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
    
//    coordinate = _bead->coordinate;
}

void Bubble::updatePositionManually() {
    //if reaching the desire position
    if(iter > SysParams::Chemistry().StepTotal) {
        iter = 1;
        currentStep++;
    }
    //All position updates will be finished in 1 second
    //Step displacement is 1 /StepTotal
    if(tau() > (currentStep * SysParams::Chemistry().StepTime + iter * 1 / SysParams::Chemistry().StepTotal)){
        double *bcoord, *coord, step;
        
        if(currentStep > SysParams::Chemistry().IterChange){
            step = SysParams::Chemistry().AFMStep2;
        }
        else{
            step = SysParams::Chemistry().AFMStep1;
        }
        coord = CUDAcommon::serlvars.coord;
        bcoord = &coord[3 * _bead->_dbIndex];
        bcoord[2] = coordinate[2] + step;
        
        vector<double> newcoord = {coordinate[0], coordinate[1], coordinate[2] + step};
        coordinate = newcoord;
        _bead->coordinate = newcoord;
        iter++;
    }

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
