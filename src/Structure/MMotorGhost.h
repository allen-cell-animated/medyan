
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------


#ifndef MEDYAN_MMotorGhost_h
#define MEDYAN_MMotorGhost_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents a cross-link between [Filaments](@ref Filament) that can move by way of
/// chemical reactions.

/*!
 *  The class describes interaction between 4 [Beads](@ref Bead) connected by a 
 *  MotorGhost, and equilibrium constants. Initial length is determinated by the condition 
 *  of zero initial stress, i.e., it is calculated within the constructor at initiation. 
 *  A ghost motor heads positions on a segment (between two consecutive beads on a filament) 
 *  determined by two numbers (0 to 1) position1 and position2 (finite number of steps 
 *  before move to the next segment o--x-o- -> o---xo- -> o---ox) held in the parent
 *  MotorGhost. They can be changed as a result of chemical reaction, and we consider that 
 *  as the motor making a step.
 */
class MMotorGhost {
    
public:
    floatingpoint stretchForce = 0.0; ///< Stretching force of motor at current state
    
    /// Main constructor
    /// @param position - position on cylinder 1 and 2, respectively
    /// @param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MMotorGhost(int motorType, int numBoundHeads, floatingpoint position1, floatingpoint position2,
                const vector<floatingpoint>& coord11, const vector<floatingpoint>& coord12,
                const vector<floatingpoint>& coord21, const vector<floatingpoint>& coord22);
    
    //@{
    /// Getter for constants
    floatingpoint getStretchingConstant(){return _kStretch;}
    floatingpoint getEqLength(){return _eqLength;}
    //@}
    
    /// Reset the spring constant of the motor based on number of bound heads
    void setStretchingConstant(int motorType, floatingpoint numBoundHeads);
    
    /// Set parent
    void setMotorGhost(MotorGhost* motor) {_pMotorGhost = motor;}
    /// Get parent
    MotorGhost* getMotorGhost() {return _pMotorGhost;}
    
    //@{
    /// Length management
    void setLength(floatingpoint l){_currentLength = l;}
    floatingpoint getLength() {return _currentLength;}
    //@}
    
private:
    floatingpoint _eqLength;  ///< Equilibrium length, set at construction
    floatingpoint _kStretch;  ///< Stretching parameter
    
    MotorGhost* _pMotorGhost; ///< Pointer to parent
    
    floatingpoint _currentLength; ///< Current length of motor
};


#endif
