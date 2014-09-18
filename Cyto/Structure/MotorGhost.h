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

#include "common.h"

class Cylinder;

//!  A MotorGhost class is used to describe mechanical properties of a motor
/*! The class describes interaction between 4 beads connected by a motor. Ghost stands for the fact that there 
 *  are NO ACTUAL BEADS assotiated with these motors, but just potentials acting on connected filament beads. Initial
 *  length of a ghost motor is determinated by the condition of zero initial stress, i.e., it calculated within the
 *  constructor at initiation. A ghost motor heads positions on a segment (between two consecutive beads on a filament) 
 *  determined by two numbers (0 to 1) position1 and position2 (finit number of steps before move to the next segment 
 *  o--x-o- -> o---xo- -> o---ox). they can be changed as a result of chemical reaction, than we consider that the motor
 *  made a step.
 */

class MotorGhost{
    
public:
    ///Main constructor
    MotorGhost(Cylinder* pc1, Cylinder* pc2, double stretchConst, double position1, double position2);
    
    //Public methods called by MotorGhost:
    double getStretching();
    
    ///Getters for mechanical properties
    Cylinder* GetFirstCylinder(){return _pc1;}
    Cylinder* GetSecondCylinder(){return _pc2;}
    double GetStretchingConstant(){return _kStretch;}
    double GetFirstPosition(){return  _position1;}
    double GetSecondPosition(){return _position2;}
    double GetEqLength(){return _eqLength;}
    
    
private:
    Cylinder* _pc1;
    Cylinder* _pc2;
    double _eqLength;
    double _kStretch;
    double _position1;
    double _position2;
};

#endif /* defined(__CytoMech__MotorGhost__) */
