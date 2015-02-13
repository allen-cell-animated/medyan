
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_RateChanger_h
#define M3SYM_RateChanger_h

#include "common.h"

/// Used to change reaction rates based on forces in the network
/*!
 *  The RateChanger class is an abstract class which allows
 *  for rate changing based on a given force. Different implementations
 *  of this class will have different rate changing models, and will
 *  all implement the changeRate() function.
 */

class RateChanger {
    
public:
    /// Change the reaction rate based on a bare rate and given force.
    virtual float changeRate(float bareRate, double force) = 0;
};

/// A brownian ratchet implementation of the DynamicRateChanger.
/// Used for filament polymerization when under load force.

/// @note - This function updates polymerization rates based on the
/// Elastic Brownian Ratchet Model (by Peskin et al, Biophys J 1993):
///
///                 k = k_0 * exp(-f * x / kT)

class BrownianRatchet : public RateChanger {
    
private:
    double _x; ///< The characteristic length for this function
public:
    BrownianRatchet(double charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A catch-slip bond implementation of the DynamicRateChanger.
///Used for either cross-linker or motor unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form (Guo et al 2006):
///
///  k = k_0 * (a_c * exp(-f * x_c / kT)  + a_s * exp(-f * x_s / kT))
///
///  where x and a are the characteristic lengths and amplitudes
///  of the catch and slip portions of the function, respectively.

class CatchSlipBond : public RateChanger {
    
private:
    double _a1; ///< catch bond amplitude
    double _a2; ///< slip bond amplitude
    double _x1; ///< catch bond characteristic length
    double _x2; ///< slip bond characteristic length
    
public:
    CatchSlipBond(double amplitude1, double amplitude2,
                  double charLength1, double charLength2) :
        _a1(amplitude1),  _a2(amplitude2),
        _x1(charLength1), _x2(charLength2) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A slip bond implementation of the DynamicRateChanger.
///Used for either cross-linker or motor unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form (Bell form, Bell et al, 1978):
///
///                 k = k_0 * exp(f * a / kT)
///
/// So as to exponetially increase the unbinding with more force.

class SlipBond : public RateChanger {
    
private:
    double _x; ///< The characteristic length for this function
    
public:
    SlipBond(double charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A stall force implementation of the DynamicRateChanger.
///Used for motor walking when under stress.

/// @note - This function updates walking rates based on the
/// following exponential form:
///
///                 k = k_0 * exp(-f * x / kT)
///
/// so as to exponentially decrease the walking rate with
/// more force.

class ExpStall : public RateChanger  {
    
private:
    double _x; ///< The characteristic length for this function
    
public:
    ExpStall(double charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

#endif
