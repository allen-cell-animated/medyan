
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

#include "SysParams.h"

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
                   double charLength1, double charLength2)
    
        : LinkerRateChanger(linkerType),
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
    BasicSlip(short linkerType, double charLength)
    
        : LinkerRateChanger(linkerType), _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A low duty PCM catch bond implementation of the MotorRateChanger
///Used for a low duty ratio motor unbinding when under stress

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates unbinding rates of a
/// Myosin II ensemble based on the following exponential form
/// (adopted from Erdmann et al 2013):
///
///      k_unbinding,eff = (k_0 / N_b) * exp(-F / (N_b * F_0))
///
/// where k_0 is the unbinding rate under zero load,
/// F_0 is the characteristic force defining this catch,
/// and N_b is the number of bound motor heads in the ensemble,
/// approximated by Erdmann et al 2013 to be:
///
///             N_b = p * N_t + (F * alpha)
///
/// where alpha has been empirically determined to be 0.04
/// for a low duty ratio motor (p = 0.1).
/// The total rate for unbinding for the entire mini-filament
/// unbinding from two actin filaments is:
///
///   k_unbinding,total = k_unbinding,eff^2 / k_binding,eff
///
/// where k_binding,eff is as previously defined.

class LowDutyPCMCatch : public MotorRateChanger {
    
private:
    double _F0;  ///< characteristic force
    
public:
    LowDutyPCMCatch(short motorType, double charForce)
    
        : MotorRateChanger(motorType), _F0(charForce) {}
    
    virtual float changeRate(float onRate, float offRate,
                             int numHeads, double force);
};


///A low duty Hill form stall force implementation of the MotorRateChanger.
///Used for a low duty ratio motor walking when under stress.

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates walking rates based on the
/// following form (based on Hill et al 1937):
///
///   k_eff = k_0 * (F_0 - F / N_t) / (F_0 + (F / (N_t * beta)))
///
/// where F_0 is the characteristic force defining this catch,
/// beta has been empirically determined to be 0.12, and
/// k_0 is the walking rate under zero load, which was approximated
/// by Erdmann et al 2013 to be:
///
///        k_0 = ((N_t - N_b) / N_b) * k_b
///
/// where k_b is the binding rate of a single motor, d_step is
/// the size of a single motor step, d_total is the total step size
/// of the ensemble in simulation, and N_t is the total number of heads.
///
/// It is noted that the true k_0 is also multipilied by a fractional
/// step size corresponding to the step size in simulation,
/// d_step / d_total where d_total is the total step size in simulation,
/// based on the number of binding sites per cylinder.
///
/// N_b is the number of bound motor heads under zero load,
/// which for the low duty ratio motor (p = 0.1) has been set to
///
///             N_b = p * Nt
///
class LowDutyHillStall : public MotorRateChanger  {
    
private:
    double _F0;  ///< characteristic force
    float _stepFrac = 1.0; ///< step size of a single head relative to sim
    
public:
    LowDutyHillStall(short motorType, short filamentType, double charForce)
    
        : MotorRateChanger(motorType), _F0(charForce) {
    
        //calculate rate based on step fraction
        double d_step = SysParams::Chemistry().motorStepSize[_motorType];
        
        double d_total = (double)SysParams::Geometry().cylinderSize[filamentType] /
                                 SysParams::Chemistry().numBindingSites[filamentType];
        
        _stepFrac = d_step / d_total;
    }
    
    virtual float changeRate(float onRate, float offRate,
                             int numHeads, double force);
};


#endif
