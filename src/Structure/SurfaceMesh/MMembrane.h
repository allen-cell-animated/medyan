#ifndef MEDYAN_MMembrane_h
#define MEDYAN_MMembrane_h

#include "SysParams.h"

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

class MMembrane {

private:

    double _eqVolume; // The equilibrium volume

    double _eqArea; // The equilibrium area, used in quadratic stretching energy
    double _kElastic; // Area elasticity, used in quadratic stretching energy

    double _tension; // Membrane tension, used in linear stretching energy.
                     // Note that the tension can also be defined in quadratic energy case,
                     // but is not stored in this variable.

public:
    MMembrane(short membraneType) {
        if(!SysParams::Mechanics().MemElasticK.empty())
            _kElastic = SysParams::Mechanics().MemElasticK[membraneType];
        if(!SysParams::Mechanics().MemTension.empty())
            _tension = SysParams::Mechanics().MemTension[membraneType];
    }

    double getEqVolume()const { return _eqVolume; }
    void setEqVolume(double eqVolume) { _eqVolume = eqVolume; }

    double getEqArea()const { return _eqArea; }
    void setEqArea(double eqArea) { _eqArea = eqArea; }

    double getKElastic()const { return _kElastic; }
    void setKElastic(double kElastic) { _kElastic = kElastic; }

    double getTension()const { return _tension; }
    void setTension(double tension) { _tension = tension; }
    
};





#endif
