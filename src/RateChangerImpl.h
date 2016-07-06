
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

#ifndef MEDYAN_RateChangerImpl_h
#define MEDYAN_RateChangerImpl_h

#include "common.h"

#include "RateChanger.h"

#include "SysParams.h"

/// A brownian ratchet implementation of the FilamentRateChanger.
/// Used for filament polymerization when under load force.

/// @note - This function updates polymerization rates based on the
/// Elastic Brownian Ratchet Model (by Peskin et al, Biophys J 1993):
///
///                 k = k_0 * exp(-f * x / kT)
///
/// @note - We note that we have implemented a load force ceiling of 100pN
///         such that recalculated reaction rates are not excessively small.
///         This would produce problems in the chemical simulation algorithm.
///         A 100pN load force ensures that the polymerization rate produced
///         will be small enough such that polymerization events are VERY rare (factor = 1E-29).

class BrownianRatchet : public FilamentRateChanger {
    
private:
    double _x; ///< The characteristic length for this function
    const double _max_f = 100; ///< 100pN ceiling
    
public:
    BrownianRatchet(double charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A catch-slip bond implementation of the LinkerRateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form of Bell et al, 1978:
///
///  k = k_0 * (a_c * exp(-f * x_c / kT)  + a_s * exp(-f * x_s / kT))
///
///  where x and a are the characteristic lengths and amplitudes
///  of the catch and slip portions of the function, respectively.

class CatchSlip : public LinkerRateChanger {
    
private:
    double _a1; ///< catch bond amplitude
    double _a2; ///< slip bond amplitude
    double _x1; ///< catch bond characteristic length
    double _x2; ///< slip bond characteristic length
    
public:
    CatchSlip(short linkerType,
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
/// following exponential form of Bell et al, 1978:
///
///                 k = k_0 * exp(f * a / kT)
///
/// So as to exponetially increase the unbinding with more force.

class Slip : public LinkerRateChanger {
    
private:
    double _x; ///< The characteristic length for this function
    
public:
    Slip(short linkerType, double charLength)
    
        : LinkerRateChanger(linkerType), _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};

///A low duty catch bond implementation of the MotorRateChanger
///Used for a low duty ratio motor unbinding when under stress
///Adopted from the Parallel Cluster Model (PCM) of Erdmann et al. 2013.

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates unbinding rates of a
/// Myosin II ensemble based on the following exponential form:
///
///      k_unbinding,eff = beta * (k_0 / N_b) * exp(-F / (N_b * F_0))
///
/// where k_0 is the unbinding rate under zero load,
/// F_0 is the characteristic force defining this catch,
/// beta has been chosen to be 0.2,
/// and N_b is the number of bound motor heads in the ensemble,
/// approximated by Erdmann et al. 2013 to be:
///
///             N_b = p * N_t + (F * gamma)
///
/// where gamma has been chosen to be 0.05
/// for a low duty ratio motor (p = 0.1).

class LowDutyCatch : public MotorRateChanger {
    
private:
    double _F0;  ///< characteristic force
    
    //@{
    ///Constant parameters
    const double dutyRatio = 0.1;
    const double beta = 0.2;
    const double gamma = 0.05;
    //@}
    
public:
    LowDutyCatch(short motorType, double charForce)
    
        : MotorRateChanger(motorType), _F0(charForce) {}
    
    virtual float changeRate(float onRate, float offRate,
                             int numHeads, double force);
};

///A low duty catch-slip bond implementation of the MotorRateChanger
///Used for a low duty ratio motor unbinding when under stress
///Adopted from both the Parallel Cluster Model (PCM) of Erdmann et al.,
///as well as Stam et al. 2015.

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates unbinding rates of a
/// Myosin II ensemble based on the following exponential form:
///
///      k_unbinding,eff = beta * (k_0 / N_b) *
///                        (_a1 * exp(-F / (N_b * _FCatch)) +
///                         _a2 * exp( F / (N_b * _FSlip))))
///
/// where k_0 is the unbinding rate under zero load,
/// _FCatch and _FSlip are the characteristic forces,
/// _a1 and _a2 are the amplitudes of each part, taken
/// as 0.92 and 0.08 respectively (Stam et al, 2015),
/// beta has been chosen to be 0.2,
/// and N_b is the number of bound motor heads in the ensemble,
/// approximated by Erdmann et al. 2013 to be:
///
///             N_b = p * N_t + (F * gamma)
///
/// where gamma has been chosen to be 0.05
/// for a low duty ratio motor (p = 0.1).

class LowDutyCatchSlip : public MotorRateChanger {
    
private:
    double _FCatch;  ///< characteristic catch force
    double _FSlip;   ///< characteristic slip force
    
    //@{
    ///Constant parameters
    const double dutyRatio = 0.1;
    const double beta = 0.2;
    const double gamma = 0.05;
    
    const double _a1 = 0.92;   ///< catch amplitude
    const double _a2 = 0.08;   ///< slip amplitude
    //@}
    
public:
    LowDutyCatchSlip(short motorType, double charCatchForce,
                                         double charSlipForce)
    
    : MotorRateChanger(motorType),
      _FCatch(charCatchForce), _FSlip(charSlipForce) {}
    
    virtual float changeRate(float onRate, float offRate,
                             int numHeads, double force);
};


///A low duty stall force implementation of the MotorRateChanger.
///Used for a low duty ratio motor walking when under stress.
///Adopted from Hill et al. 1937, and Erdmann et al. 2013.

/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates walking rates based on the Hill form:
///
///   k_eff = k_0 * (F_0 - F / N_t) / (F_0 + (F / (N_t * zeta)))
///
/// where F_0 is the characteristic force defining this stall,
/// zeta has been chosen to be 0.1, and
/// k_0 is the walking rate under zero load, which was approximated
/// by Erdmann et al. 2013 to be:
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
class LowDutyStall : public MotorRateChanger  {
    
private:
    double _F0;            ///< characteristic force
    float _stepFrac = 1.0; ///< step size of a single head relative to sim
    
    
    //@{
    ///Constant parameters
    const double dutyRatio = 0.1;
    const double zeta = 0.1;
    //@}
    
    
public:
    LowDutyStall(short motorType, short filamentType, double charForce)
    
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
