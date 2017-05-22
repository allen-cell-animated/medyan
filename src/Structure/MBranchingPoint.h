
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

#ifndef MEDYAN_MBranchingPoint_h
#define MEDYAN_MBranchingPoint_h

#include "common.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents the mechanical component of a BranchingPoint.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a 
 *  BranchingPoint, and its associated equilibrium constants. The branch is a connection 
 *  between two [Cylinders](@ref Cylinder), where one resides on the mother Filament and
 *  the other is the actual branch. This is positioned on a Cylinder with the field
 *  position (between 0 and 1), held by the parent BranchingPoint.
 */
class MBranchingPoint {
    
public:
    /// Main constructor, sets constants
    MBranchingPoint(int branchType);
    
    //@{
    /// Getter for constants
    double getStretchingConstant(){return _kStretch;}
    double getEqLength(){return _eqLength;}
    
    double getBendingConstant(){return _kBend;}
    double getEqTheta(){return _eqTheta;}
    
    double getDihedralConstant(){return _kDihedr;}
    
    double getPositionConstant(){return _kPosition;}
    //@}
    
    /// Set parent
    void setBranchingPoint(BranchingPoint* BranchingPoint) {
        _pBranchingPoint = BranchingPoint;
    }
    /// Get parent
    BranchingPoint* getBranchingPoint() {return _pBranchingPoint;}
    
    //Qin -----
    double stretchForce = 0.0; ///< Stretching force of brancher at current state
    
private:
    double _eqLength;  ///< Equilibrium length
    double _kStretch;  ///< Stretching constant
    
    double _eqTheta; ///< Bending equilibrium angle
    double _kBend; //< Bending constant
    
    double _kDihedr; ///< Twisting constant
    
    double _kPosition; ///< Position constant
    
    BranchingPoint* _pBranchingPoint; ///< Pointer to parent branch point
    
};


#endif
