//
//  MCylinder.h
//  CytoMech
//
//  Created by Konstantin Popov on 6/30/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MCylinder__
#define __CytoMech__MCylinder__

#include <iostream>
#include <vector>

#include "common.h"

class Cylinder;
class Filament;
class Bead;

///MCylinder class is used to hold mechanical properties of a cylinder
/*!
 * MCylinder is a class to hold mechanical properties of a cylinder, including equilibrium force constants and
 * pointers to related structures, including a neighbors list of cylinders as well as its corresponding beads.
 */

class MCylinder {

private:
    Bead* _pFirst;  ///< Pointer to the first bead, associated with this cylinder ;
    Bead* _pSecond; ///< Pointer to the end bead in the cylinder. Either empty - last cylinder, or pointer to the first Bead in a next cylinder.
    
    Cylinder* _pCylinder;
    
// Mechanical constants:
    
    double _eqLength;   //< Lenght of unstertched cylinder;
    double _eqAngle;   //< Equilibrium value for angle in bending potential. For interaction between this cylinder and PREVIOUS;
    double _eqAngleTwist;
    double _kStretch;  //< Local stretching constant, describes axial stretching of a single cylinder;
    double _kBend;  //< Local bending constant, which describes bendinig interaction between current and PREVIOUS cylinders;
    double _kTwist;
    
    std::vector<MCylinder*> _NeighbourList;

public:
    
    ///Constructor and destructor
    MCylinder(Filament* pf, Bead* firstBead, Bead* secondBead, double eqLength);
    virtual ~MCylinder() {}
    
    ///Other setter and getter functions:
    
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    Cylinder* getCylinder() {return _pCylinder;}
    
    Bead* GetFirstBead();
    Bead* GetSecondBead();
    
    void SetEqLength(double L);
    double GetEqLength();
    
    void SetAngle(double alpha);
    double GetAngle();
    
    void SetStretchingConst(double k);
    double GetStretchingConst();
    
    void SetBendingConst(double k);
    double GetBendingConst();
    
    void SetTwistingConst(double k);
    double GetTwistingConst();
    
//  void DeleteBead(Bead*) {std::cout<<"not implemented"<<std::endl;}
    
    
};

#endif /* defined(__CytoMech__MCylinder__) */
