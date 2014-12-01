//
//  MMotorGhost.h
//  Cyto
//
//  Created by James Komianos on 10/20/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MMotorGhost__
#define __Cyto__MMotorGhost__

#include <iostream>
#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

///MMotorGhost class represents a cross-link between filaments that can move across a filament

 /*! The class describes interaction between 4 beads connected by a motor. Ghost stands for the fact that there
 *  are NO ACTUAL BEADS assotiated with these motors, but just potentials acting on connected filament beads. Initial
 *  length of a ghost motor is determinated by the condition of zero initial stress, i.e., it calculated within the
 *  constructor at initiation. A ghost motor heads positions on a segment (between two consecutive beads on a filament)
 *  determined by two numbers (0 to 1) position1 and position2 (finit number of steps before move to the next segment
 *  o--x-o- -> o---xo- -> o---ox). they can be changed as a result of chemical reaction, than we consider that the motor
 *  made a step.
 */

class MMotorGhost {
    
public:
    ///Main constructor
    ///@param position - position on cylinder 1 and 2, respectively
    ///@param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MMotorGhost(double stretchConst, double position1, double position2,
            const vector<double>& coord11, const vector<double>& coord12,
            const vector<double>& coord21, const vector<double>& coord22);
    
    ///Getters for constants
    double getStretchingConstant(){return _kStretch;}
    double getEqLength(){return _eqLength;}
    
    ///Setter and getter for parent linker
    void setMotorGhost(MotorGhost* motor) {_pMotorGhost = motor;}
    MotorGhost* getMotorGhost() {return _pMotorGhost;}
    
private:
    double _eqLength;
    double _kStretch;
    
    MotorGhost* _pMotorGhost; ///< ptr to parent linker
    
};


#endif /* defined(__Cyto__MMotorGhost__) */
