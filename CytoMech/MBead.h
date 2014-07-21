//
//  MBead.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MBead__
#define __CytoMech__MBead__

#include <iostream>
#include <vector>
#include <list>
#include "Mcommon.h"
#include "MComponent.h"


class Bead : public MComponent
{
public:
    
    ///The basic and simplest mecanical class. Contains information about a bead and some mecanical constant needed to calculate interactions. Constants are local and related to a bond between i(curent) and i-1 bead.
    
    Bead (std::vector<double> v): coordinate(v), _parent(NULL), force(3, 0), forceAux(3, 0) {}
    Bead(): coordinate (3, 0), _parent(NULL), force(3, 0), forceAux(3, 0) {}
//    virtual ~Bead();
    
	std::vector<double> coordinate; //Coordinates of the bead
	std::vector<double> force;      // Forces based on curent coordinates. Forces should alwais corespond to current coordinates.
    std::vector<double> forceAux; //An Aux field neede during CG minimization.
    
    
    //Aux functios:
    double CalcForceSquare() {return force[0]*force[0] + force[1]*force[1] + force[2]*force[2];  } // Aux method for CG minimization;
    double CalcForceSquare( int i) {return forceAux[0]*forceAux[0] + forceAux[1]*forceAux[1] + forceAux[2]*forceAux[2];  }
    double CalcDotForceProduct() { return force[0]*forceAux[0] + force[1]*forceAux[1] + force[2]*forceAux[2];}
    
    //Function which returns a pointer to a filament(parent "Composite pattern") to which this bead belongs:
    MComposite* GetParent(){
        return _parent;
    }
    
    void SetParent(MComposite* pc) {_parent = pc;}
    
private:
    
    MComposite *_parent; /// A pointer to a parent class (Filament, linker, motor)
};


#endif /* defined(__CytoMech__MBead__) */
