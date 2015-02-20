
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

/// A brownian ratchet implementation of the RateChanger.
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

///A catch-slip bond implementation of the RateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form (Guo et al 2006):
///
///  k = k_0 * (a_c * exp(-f * x_c / kT)  + a_s * exp(-f * x_s / kT))
///
///  where x and a are the characteristic lengths and amplitudes
///  of the catch and slip portions of the function, respectively.

class LinkerCatchSlip : public RateChanger {
    
private:
    double _a1; ///< catch bond amplitude
    double _a2; ///< slip bond amplitude
    double _x1; ///< catch bond characteristic length
    double _x2; ///< slip bond characteristic length
    
public:
    LinkerCatchSlip(double amplitude1, double amplitude2,
                    double charLength1, double charLength2) :
        _a1(amplitude1),  _a2(amplitude2),
        _x1(charLength1), _x2(charLength2) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A slip bond implementation of the RateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form (Bell form, Bell et al, 1978):
///
///                 k = k_0 * exp(f * a / kT)
///
/// So as to exponetially increase the unbinding with more force.

class LinkerSlip : public RateChanger {
    
private:
    double _x; ///< The characteristic length for this function
    
public:
    LinkerSlip(double charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A PCM catch bond implementation of the RateChanger
///Used for motor unbinding when under stress

/// @note - This function updates unbinding rates based on the
/// following exponential form (Erdmann et al, JACS 2013):
///
///      k_eff = N_b * k_unbinding * exp(-F / (N_b * F_0))
///
/// where F_0 is the characteristic force defining this catch,
/// and N_b is the number of bound motor heads determined by
///
///             N_b = 0.33*N_t + (F * alpha)
///
/// where alpha has been empirically determined to be 0.028
class MotorPCMCatch : public RateChanger {
    
private:
    
    double _F0;  ///< characteristic force
    int _Nt; ///< number of motor heads
    
public:
    MotorPCMCatch(double charForce, int numHeads)
        : _F0(charForce), _Nt(numHeads) {}
    
    virtual float changeRate(float bareRate, double force);
};


///A Hill form stall force implementation of the RateChanger.
///Used for motor walking when under stress.

/// @note - This function updates walking rates based on the
/// following form (based on Erdmann et al, JACS 2013):
///
///   k_eff = (k_w / N_b) * (F_0 - F) / (F_0 + (F / beta))
///
/// where F_0 is the characteristic force defining this catch,
/// beta has been empirically determined to be 0.21, and
/// N_b is the number of bound motor heads determined by
///
///             N_b = 0.33*N_t + (F * alpha)
///
/// where alpha has been empirically determined to be 0.028

class MotorHillStall : public RateChanger  {
    
private:
    double _F0;  ///< characteristic force
    int _Nt; ///< number of motor heads
    
public:
    MotorHillStall(double charForce, int numHeads)
        : _F0(charForce), _Nt(numHeads) {}
    
    virtual float changeRate(float bareRate, double force);
};

#endif
