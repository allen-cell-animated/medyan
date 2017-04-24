
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

class LinkerCatchSlip : public LinkerRateChanger {
    
private:
    double _a1; ///< catch bond amplitude
    double _a2; ///< slip bond amplitude
    double _x1; ///< catch bond characteristic length
    double _x2; ///< slip bond characteristic length
    
public:
    LinkerCatchSlip(short linkerType,
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

class LinkerSlip : public LinkerRateChanger {
    
private:
    double _x; ///< The characteristic length for this function
    
public:
    LinkerSlip(short linkerType, double charLength)
    
    : LinkerRateChanger(linkerType), _x(charLength) {}
    
    virtual float changeRate(float bareRate, double force);
};


///A catch bond implementation of the MotorRateChanger
///Used for a motor unbinding when under stress
///Adopted from the results of Erdmann et al. 2013.

/// @note - This function updates unbinding rates of a
/// Myosin II ensemble based on the following exponential form:
///
///    k_unbinding,eff = k_0 * exp(-F / F_0)
///
/// where k_0 is the unbinding rate under zero load, F_0 is the characteristic
/// force defining this catch.

class MotorCatch : public MotorRateChanger {
    
private:
    double _k_0 = 0.2; ///< unbinding rate bare (/s)
    double _F0_slip = 100; ///< slip force
    
public:
    double _F0_catch = 25;  ///< catch force
    
    MotorCatch(short motorType)
    
    : MotorRateChanger(motorType) {}
    
    virtual float changeRate(float onRate, float offRate, double force);
};

///Implementation1
class MotorACatch : public MotorCatch {
    
public:
    MotorACatch(short motorType)
    
    : MotorCatch(motorType){
        _F0_catch = _F0_catch * SysParams::Chemistry().sigma;
    }
    
    
};

///Implementation2
class MotorBCatch : public MotorCatch {
    
public:
    MotorBCatch(short motorType)
    
    : MotorCatch(motorType){}
};


///A stall force implementation of the MotorRateChanger.
///Used for a motor walking when under stress.
///Adopted from Hill et al. 1937

/// @note - This function updates walking rates based on a simplified Hill form:
///
///   k_eff = k_0 * (1 - F/F_0)
///
/// where F_0 is the characteristic force defining this stall,
/// k_0 is the walking rate under zero load
///
class MotorStall : public MotorRateChanger  {
    
private:
    double _F0 = 20;       ///< characteristic force in pN
    double _stepFrac = 1.0; ///< step size of a single head relative to sim
    
public:
    //FOR MYOSIN-ISOFORMS
    double v_0 = 20;
    
    MotorStall(short motorType, short filamentType)
    
    : MotorRateChanger(motorType) {
        
        double d_total = (double)SysParams::Geometry().cylinderSize[filamentType] /
        SysParams::Chemistry().numBindingSites[filamentType];
        
        _stepFrac = d_total;
    }
    
    //virtual float numBoundHeads(float onRate, float offRate,
    //                            double force, int numHeads) {return 0.0;}
    
    virtual float changeRate(float onRate, float offRate, double force);
};

///Implementation1
class MotorAStall : public MotorStall {
    
public:
    MotorAStall(short motorType, short filamentType)
    
    : MotorStall(motorType, filamentType){
        v_0 = v_0 * SysParams::Chemistry().lambda;
    }
};

///Implementation2
class MotorBStall : public MotorStall {
    
public:
    MotorBStall(short motorType, short filamentType)
    
    : MotorStall(motorType, filamentType){}
};



#endif
