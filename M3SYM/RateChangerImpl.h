
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

#ifndef M3SYM_RateChangerImpl_h
#define M3SYM_RateChangerImpl_h

#include "common.h"

#include "RateChanger.h"

/// A brownian ratchet implementation of the FilamentRateChanger.
/// Used for filament polymerization when under load force.

/// @note - This function updates polymerization rates based on the
/// Elastic Brownian Ratchet Model (by Peskin et al, Biophys J 1993):
///
///                 k = k_0 * exp(-f * x / kT)

class BrownianRatchet : public FilamentRateChanger {
    
private:
    double _x; ///< The characteristic length for this function
public:
    BrownianRatchet(double charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A catch-slip bond implementation of the LinkerRateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form (Bell form, Bell et al, 1978):
///
///  k = k_0 * (a_c * exp(-f * x_c / kT)  + a_s * exp(-f * x_s / kT))
///
///  where x and a are the characteristic lengths and amplitudes
///  of the catch and slip portions of the function, respectively.

class BasicCatchSlip : public LinkerRateChanger {
    
private:
    double _a1; ///< catch bond amplitude
    double _a2; ///< slip bond amplitude
    double _x1; ///< catch bond characteristic length
    double _x2; ///< slip bond characteristic length
    
public:
    BasicCatchSlip(short linkerType,
                   double amplitude1, double amplitude2,
                   double charLength1, double charLength2) :
    LinkerRateChanger(linkerType),
    _a1(amplitude1),  _a2(amplitude2),
    _x1(charLength1), _x2(charLength2) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A slip bond implementation of the LinkerRateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form (Bell form, Bell et al, 1978):
///
///                 k = k_0 * exp(f * a / kT)
///
/// So as to exponetially increase the unbinding with more force.

class BasicSlip : public LinkerRateChanger {
    
private:
    double _x; ///< The characteristic length for this function
    
public:
    BasicSlip(short linkerType, double charLength) :
    LinkerRateChanger(linkerType), _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A low duty PCM catch bond implementation of the MotorRateChanger
///Used for a low duty ratio motor unbinding when under stress

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates unbinding rates based on the
/// following exponential form (Erdmann et al, JACS 2013):
///
///      k_eff = (k_0 / N_b) * exp(-F / (N_b * F_0))
///
/// where k_0 is the unbinding rate under zero load,
/// F_0 is the characteristic force defining this catch,
/// and N_b is the number of bound motor heads in the ensemble,
/// approximated by Erdmann et al 2013 to be:
///
///             N_b = 0.1 * N_t + (F * alpha)
///
/// where alpha has been empirically determined to be 0.04
/// for a low duty ratio motor (p = 0.1).

class LowDutyPCMCatch : public MotorRateChanger {
    
private:
    double _F0;  ///< characteristic force
    
public:
    LowDutyPCMCatch(short motorType, double charForce) :
    MotorRateChanger(motorType), _F0(charForce) {}
    
    virtual float changeRate(float bareRate, int numHeads, double force);
};


///A low duty Hill form stall force implementation of the MotorRateChanger.
///Used for a low duty ratio motor walking when under stress.

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates walking rates based on the
/// following form (based on Hill et al, Erdmann et al 2013):
///
///   k_eff = k_0 * (F_0 - F / N_t) / (F_0 + (F / (N_t * beta)))
///
/// where F_0 is the characteristic force defining this catch,
/// beta has been empirically determined to be 0.12, and
/// k_0 is the walking rate under zero load, which was approximated
/// by Erdmann et al 2013 to be:
///
///  k_0 = ((N_t - N_b) / N_b) * (d_step / d_total) * k_b
///
/// where k_b is the binding rate of a single motor, d_step is
/// the size of a single motor step, d_total is the total step size
/// of the ensemble in simulation, and N_t is the total number of heads.
///
/// N_b is the number of bound motor heads under zero load,
/// which for the low duty ratio motor (p = 0.1) has been set to
///
///             N_b = 0.1 * Nt
///
class LowDutyHillStall : public MotorRateChanger  {
    
private:
    double _F0;  ///< characteristic force
    
public:
    LowDutyHillStall(short motorType, double charForce) :
    MotorRateChanger(motorType), _F0(charForce) {}
    
    virtual float changeRate(float bareRate, int numHeads, double force);
};


#endif