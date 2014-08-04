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

class MCylinder
{

private:
    Bead* _pFirst;  //< Pointer to the first bead, associated with this cylinder ;
    Bead* _pSecond; //< Pointer to the end bead in the cylinder. Either empty - last cylinder, or pointer to the first Bead in a next cylinder.
    
    
// Mechanical constants:
    double _eqLenght;   //< Lenght of unstertched cylinder;
    double _eqAngle;   //< Equilibrium value for angale in bending potential. For interaction between this cylinder and PREVIOUS;
    double _eqAngleTwist;
    double _kStretch;  //< Local stretching constant, describes axial stretching of a single cylinder;
    double _kBend;  //< Local bending constant, which describes bendinig interaction between current and PREVIOUS cylinders;
    double _kTwist;
    


    std::vector<MCylinder*> _NeighbourList;

public:
    
//    Cylinder();
//    Cylinder(Filament* pf, Bead* firstBead, double L, double theta, double ks, double kb, double kt);
    MCylinder(Filament* pf, Bead* firstBead);
    virtual ~MCylinder() {}
    
   
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
    
//  void DeleteBead(Bead*) {std::cout<<"not implemented"<<std::endl;}
    
    
};

#endif /* defined(__CytoMech__MCylinder__) */
