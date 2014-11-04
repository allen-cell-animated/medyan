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
    Cylinder* _pCylinder;
    
/// Mechanical constants:
    double _eqLength;   ///< Lenght of unstertched cylinder;
    double _eqAngle;   ///< Equilibrium value for angle in bending potential. For interaction between this cylinder and PREVIOUS;
    double _kStretch;  ///< Local stretching constant, describes axial stretching of a single cylinder;
    double _kBend;  ///< Local bending constant, which describes bending interaction between current and PREVIOUS cylinders;
    double _kTwist; ///< Local twisting constant, which describes stretching interaction between current and PREVIOUS cylinders
    double _kExVol; ///< Local excluded volume constant, which describes excluded volume interactions between cylinders
    
    std::vector<MCylinder*> _exVolNeighborsList; ///< list of interacting cylinders for excluded volume.

    std::vector<double> _coordinate; ///< coordinate of this MCylinder
public:
    
    ///Constructor and destructor
    MCylinder(double eqLength);
    ~MCylinder() {}
    
    ///update excluded volume neighbors list
    void updateExVolNeighborsList(std::vector<MCylinder*>& nearbyMCylinders);
    
    ///Other setter and getter functions:
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    Cylinder* getCylinder() {return _pCylinder;}
    
    void setCoordinate(std::vector<double> coordinate) {_coordinate = coordinate;}
    
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
    
    void SetExVolConst(double k);
    double GetExVolConst();
    
    std::vector<MCylinder*>& getExVolNeighborsList() {return _exVolNeighborsList;}
    
};

#endif /* defined(__CytoMech__MCylinder__) */
