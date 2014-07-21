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
#include "Mcommon.h"
#include "MComposite.h"

class Cylinder : public MComposite
{

private:
    Filament* _pFilament; //< Pointer to filament where this cilinder belongs;
    Bead* _pFirst;  //< Pointer to the first bead, associated with this cylinder ;
    Bead* _pSecond; //< Pointer to the end bead in the cylinder. Either empty - last cylinder, or pointer to the first Bead in a next cylinder.
    
    int _positionFilament; // position on filament (1st, 2nd, ... etc);
    
    
// Mechanical constants:
    double _eqLenght;   //< Lenght of unstertched cylinder;
    double _eqAngle;   //< Equilibrium value for angale in bending potential. For interaction between this cylinder and PREVIOUS;
    double _eqAngleTwist;
    double _kStretch;  //< Local stretching constatn, describes axial stretching of a single cylinder;
    double _kBend;  //< Local bending constant, which describes bendinig interaction between current and PREVIOUS cylinders;
    double _kTwist;
    
    bool _ifLast;

    std::vector<Cylinder*> _NeighbourList;

public:
    
//    Cylinder();
//    Cylinder(Filament* pf, Bead* firstBead, double L, double theta, double ks, double kb, double kt);
    Cylinder(Filament* pf, Bead* firstBead);
    virtual ~Cylinder() {}
    
    bool IfLast();
    void SetLast(bool);
    void SetSecondBead(Bead *pb);
    
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
    
    void DeleteBead(Bead*) {std::cout<<"not implemented"<<std::endl;}
    
    
};

#endif /* defined(__CytoMech__MCylinder__) */
