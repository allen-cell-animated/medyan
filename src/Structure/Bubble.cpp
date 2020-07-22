
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
#include "Structure/Special/AFM.h"

#include "SysParams.h"
#include "CUDAcommon.h"

Bubble::Bubble(SubSystem* ps, vector<floatingpoint> coordinates, short type)

    : Trackable(true, false, true, false), _ps(ps), _type(type),
      coord(vector2Vec< 3, floatingpoint >(coordinates)), force {},
      birthTime_(tau())
{
    
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
          
      if(SysParams::Mechanics().AFMBendingK.size() == 0){
          _AFMBendingK = 0.0;
      }
      else{
          _AFMBendingK = SysParams::Mechanics().AFMBendingK[_type];
      }
          
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
        floatingpoint step;
        
        if(currentStep > SysParams::Chemistry().IterChange){
            step = SysParams::Chemistry().AFMStep2;
        }
        else{
            step = SysParams::Chemistry().AFMStep1;
        }

        coord[2] += step;

        // Update boundary element coordinate
        static_cast<AFM*>(getParent())->getPlaneBoundaryElement()->updateCoords(vec2Vector(coord));

        iter++;
    }

}

void Bubble::printSelf()const {
    
    cout << endl;
    
    cout << "Bubble: ptr = " << this << endl;
    cout << "Bubble ID = " << getId() << endl;
    cout << "Bubble type = " << _type << endl;
    cout << "Bubble coord = " << coord << endl;
    cout << "Bubble force = " << force << endl;
    cout << "Bubble radius = " << _radius << endl;
    
    cout << endl;
}
