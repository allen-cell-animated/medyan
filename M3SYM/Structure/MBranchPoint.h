
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_MBranchPoint_h
#define M3SYM_MBranchPoint_h

#include "common.h"

//FORWARD DECLARATIONS
class BranchPoint;

/// Represents the mechanical component of a BranchPoint.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a BranchPoint, and its 
 *  associated equilibrium constants. The branch is a connection between two [Cylinders](@ref cylinder), 
 *  where one resides on the mother Filament and the other is the actual branch.
 */
class MBranchPoint {
    
public:
    /// Main constructor, sets constants
    MBranchPoint(double stretchConst, double eqLength,
                 double bendConst, double bendTheta,
                 double twistConst, double twistAngle)
                 : _kStretch(stretchConst), _eqLength(eqLength), _kBend(bendConst),
                   _eqTheta(bendTheta), _kTwist(twistConst), _eqPhi(twistAngle) {}
    
    //@{
    /// Getter for constants
    double getStretchingConstant(){return _kStretch;}
    double getEqLength(){return _eqLength;}
    
    double getBendingConstant(){return _kBend;}
    double getEqTheta(){return _eqTheta;}
    
    double getTwistingConstant(){return _kTwist;}
    double getEqPhi(){return _eqPhi;}
    //@}
    
    /// Set parent
    void setBranchPoint(BranchPoint* branchPoint) {_pBranchPoint = branchPoint;}
    /// Get parent
    BranchPoint* getBranchPoint() {return _pBranchPoint;}
    
private:
    double _eqLength;  ///< Equilibrium length
    double _kStretch;  ///< Stretching constant
    
    double _eqTheta; ///< Bending equilibrium angle
    double _kBend; //< Bending constant
    
    double _eqPhi; ///< Twisiting angle
    double _kTwist; ///< Twisting constant
    
    BranchPoint* _pBranchPoint; ///< Pointer to parent branch point
    
};


#endif
