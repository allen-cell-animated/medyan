
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

#ifndef M3SYM_MCylinder_h
#define M3SYM_MCylinder_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Cylinder;

///MCylinder class is used to hold mechanical properties of a [Cylinder](@ref Cylinder)
/*!
 * MCylinder is a class to hold mechanical properties of a [Cylinder](@ref Cylinder), including equilibrium force constants.
 */
class MCylinder {

private:
    Cylinder* _pCylinder;  ///< parent cylinder

    double _eqLength; ///< Length of unstretched cylinder
    double _eqAngle;  ///< Equilibrium value for angle in bending potential. For interaction between this cylinder and PREVIOUS
    double _kStretch; ///< Local stretching constant, describes axial stretching of a single cylinder
    double _kBend;  ///< Local bending constant, which describes bending interaction between current and PREVIOUS cylinders
    double _kTwist; ///< Local twisting constant, which describes stretching interaction between current and PREVIOUS cylinders
    double _kExVol; ///< Local excluded volume constant, which describes excluded volume interactions between cylinders
    
public:
    MCylinder(double eqLength);
    ~MCylinder() {};

    /// Set parent cylinder
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    Cylinder* getCylinder() {return _pCylinder;}
    
    /// Set the equlilibrium length, which changes mechanical constants accordingly
    void setEqLength(double L);
    /// Get the current equlibrium length of this MCylinder
    double getEqLength() {return _eqLength;}
    
    //@{
    /// Mechanical parameter management function
    void setAngle(double alpha) {_eqAngle = alpha;}
    double getAngle() {return _eqAngle;}
    
    void setStretchingConst(double k) {_kStretch = k;}
    double getStretchingConst() {return _kStretch;}
    
    void setBendingConst(double k) {_kBend = k;}
    double getBendingConst() {return _kBend;}
    
    void setTwistingConst(double k) {_kTwist = k;}
    double getTwistingConst() {return _kTwist;}
    
    void setExVolConst(double k) {_kExVol = k;}
    double getExVolConst() {return _kExVol;}
    //@}
    
};

#endif
