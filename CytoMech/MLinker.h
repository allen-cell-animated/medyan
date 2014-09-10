//
//  MLinker.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MLinker__
#define __CytoMech__MLinker__

#include <iostream>
#include "Cylinder.h"


///Linker class represents a cross-link between filaments
/*!
 *  A linker represents a cross link between cylinders. It contains mechanical information, including
 *  constants and pointers to its respective cylinders.
 */
class Linker {
    
public:
    ///Main constructor
    Linker(Cylinder* pc1, Cylinder* pc2, double stretchConst, double position1, double position2);
    
    ///Getters for constants and cylinders
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

#endif /* defined(__CytoMech__MLinker__) */
