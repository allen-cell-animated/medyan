
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
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
    double stretchForce = 0.0; ///< Stretching force of motor at current state
    
    /// Main constructor
    /// @param position - position on cylinder 1 and 2, respectively
    /// @param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MMotorGhost(int motorType, int numHeads, double position1, double position2,
                const vector<double>& coord11, const vector<double>& coord12,
                const vector<double>& coord21, const vector<double>& coord22);
    
    //@{
    /// Getter for constants
    double getStretchingConstant(){return _kStretch;}
    double getEqLength(){return _eqLength;}
    //@}
    
    /// Set parent
    void setMotorGhost(MotorGhost* motor) {_pMotorGhost = motor;}
    /// Get parent
    MotorGhost* getMotorGhost() {return _pMotorGhost;}
    
    //@{
    /// Length management
    void setLength(double l){_currentLength = l;}
    double getLength() {return _currentLength;}
    //@}
    
private:
    double _eqLength;  ///< Equilibrium length, set at construction
    double _kStretch;  ///< Stretching parameter
    
    MotorGhost* _pMotorGhost; ///< Pointer to parent
    
    double _currentLength; ///< Current length of motor
};


#endif
