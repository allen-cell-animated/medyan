//
//  MotorGhost.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MotorGhost__
#define __CytoMech__MotorGhost__

#include <iostream>
#include "Composite.h"
#include "CMotorGhost.h"
#include "MMotorGhost.h"
#include "common.h"

///FORWARD DECLARATIONS
class Cylinder;

///MotorGhost class is a wrapper for a chemical and mechanical motor
/*!
 * MotorGhost class is used to create a chemical and mechanical motor when needed.
 * It contains a constructor as well as getters for mmotorghost and cmotorghost.
 */

class MotorGhost : public Composite{
   
private:
    unique_ptr<MMotorGhost> _mMotorGhost; ///< ptr to mMotorGhost
    unique_ptr<CMotorGhost> _cMotorGhost; ///< ptr to cMotorGhost
    
    Cylinder* _c1; ///< first cylinder the linker is bound to
    Cylinder* _c2; ///< second cylinder the linker is bound to
    
    double _position1; ///< position on first cylinder
    double _position2; ///< position on second cylinder
    
    short _motorType; ///integer specifying the type of linker
    int _motorID; ///integer ID of this motor
    
    float _birthTime; ///< birth time of this motor
    
    Compartment* _compartment; ///< Compartment that this linker is in
    
    
public:
    vector<double> coordinate; ///< coordinate of midpoint, updated with updatePosition()
    
    MotorGhost(Cylinder* c1, Cylinder* c2, short motorType, int motorID, double position1, double position2, bool creation);
    ~MotorGhost();
    
    ///get cylinders
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    
    ///set cylinders
    void setFirstCylinder(Cylinder* cylinder) {_c1 = cylinder;}
    void setSecondCylinder(Cylinder* cylinder) {_c2 = cylinder;}
    
    ///setter for mlinkers and clinkers
    void setCMotorGhost(CMotorGhost* cMotorGhost) {_cMotorGhost = unique_ptr<CMotorGhost>(cMotorGhost);}
    CMotorGhost* getCMotorGhost() {return _cMotorGhost.get();}
    
    MMotorGhost* getMMotorGhost() {return _mMotorGhost.get();}
    
    ///Getters and setters for position
    double getFirstPosition() {return _position1;}
    void setFirstPosition(double position1) {_position1 = position1;}
    double getSecondPosition() {return _position2;}
    void setSecondPosition(double position2) {_position2 = position2;}
    
    short getMotorType() {return _motorType;}
    int getMotorID() {return _motorID;}
    
    ///Update the position of this Linker
    ///@note - changes compartment of clinker if needed
    void updatePosition();


};

#endif /* defined(__CytoMech__MotorGhost__) */
