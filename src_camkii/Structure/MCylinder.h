
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

#ifndef MEDYAN_MCylinder_h
#define MEDYAN_MCylinder_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Cylinder;

/// Used to hold mechanical properties of a Cylinder.

/*! The class describes interaction between 2 [Beads](@ref Bead) connected by a
 *  Cylinder in the same Filament, and its associated equilibrium constants. 
 */
class MCylinder {

private:
    Cylinder* _pCylinder;  ///< parent cylinder

    double _eqLength; ///< Length of unstretched cylinder
    double _eqTheta;  ///< Equilibrium value for angle in bending potential.
                      ///< For interaction between this cylinder and PREVIOUS
    double _eqPhi;    ///< Equilibrium value of twisiting potential
    double _kStretch; ///< Local stretching constant, describes axial stretching
                      ///< of a single cylinder
    double _kBend;    ///< Local bending constant, which describes bending
                      ///< interaction between current and PREVIOUS cylinders
    double _kTwist;   ///< Local twisting constant, which describes stretching
                      ///< interaction between current and PREVIOUS cylinders
    double _kExVol;   ///< Local excluded volume constant, which describes
                      ///< excluded volume interactions between cylinders
    
    double _currentLength; ///< The current length of the cylinder
    
public:
    /// Constructor sets equlilibrium length, and also adjusts other
    /// parameters according to this length
    MCylinder(short filamentType, double eqLength);
    ~MCylinder() {};

    /// Set parent 
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    Cylinder* getCylinder() {return _pCylinder;}
    
    /// Set the equlilibrium length, which changes mechanical constants accordingly
    void setEqLength(short filamentType, double l);
    /// Get the current equlibrium length of this MCylinder
    double getEqLength() {return _eqLength;}
    
    //@{
    /// Mechanical parameter management function
    void setEqTheta(double theta) {_eqTheta = theta;}
    double getEqTheta() {return _eqTheta;}
    
    void setEqPhi(double phi) {_eqPhi = phi;}
    double getEqPhi() {return _eqPhi;}
    
    void setStretchingConst(double k) {_kStretch = k;}
    double getStretchingConst() {return _kStretch;}
    
    void setBendingConst(double k) {_kBend = k;}
    double getBendingConst() {return _kBend;}
    
    void setTwistingConst(double k) {_kTwist = k;}
    double getTwistingConst() {return _kTwist;}
    
    void setExVolConst(double k) {_kExVol = k;}
    double getExVolConst() {return _kExVol;}
    
    void setLength(double l){_currentLength = l;}
    double getLength() {return _currentLength;}
    //@}
    
};

#endif