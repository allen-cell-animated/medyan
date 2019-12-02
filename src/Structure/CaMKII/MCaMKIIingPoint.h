
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

#ifndef MEDYAN_MCaMKIIingPoint_h
#define MEDYAN_MCaMKIIingPoint_h

#include "common.h"

//FORWARD DECLARATIONS
class CaMKIIingPoint;

/// Represents the mechanical component of a CaMKIIingPoint.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a 
 *  CaMKIIingPoint, and its associated equilibrium constants. The camkii is a connection 
 *  between two [Cylinders](@ref Cylinder), where one resides on the mother Filament and
 *  the other is the actual camkii. This is positioned on a Cylinder with the field
 *  position (between 0 and 1), held by the parent CaMKIIingPoint.
 */
class MCaMKIIingPoint {
    
public:
    /// Main constructor, sets constants
    MCaMKIIingPoint(int camkiiType);
    
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
    void setCaMKIIingPoint(CaMKIIingPoint* CaMKIIingPoint) {_pCaMKIIingPoint = CaMKIIingPoint;}
    /// Get parent
    CaMKIIingPoint* getCaMKIIingPoint() {return _pCaMKIIingPoint;}
    
private:
    double _eqLength;  ///< Equilibrium length
    double _kStretch;  ///< Stretching constant
    
    double _eqTheta; ///< Bending equilibrium angle
    double _kBend; //< Bending constant
    
    double _kDihedr; ///< Twisting constant
    
    double _kPosition; ///< Position constant
    
    CaMKIIingPoint* _pCaMKIIingPoint; ///< Pointer to parent camkii point
    
};


#endif
