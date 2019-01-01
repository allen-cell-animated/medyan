#ifndef MEDYAN_MMembrane_h
#define MEDYAN_MMembrane_h

#include "SysParams.h"

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

class MMembrane {

private:

    double _eqVolume; // The equilibrium volume

    double _eqArea; // The equilibrium area
    double _kElastic; // Area elasticity

public:
    MMembrane(short membraneType) {
        if(!SysParams::Mechanics().MemElasticK.empty())
            _kElastic = SysParams::Mechanics().MemElasticK[membraneType];
    }

    double getEqVolume()const { return _eqVolume; }
    void setEqVolume(double eqVolume) { _eqVolume = eqVolume; }

    double getEqArea()const { return _eqArea; }
    void setEqArea(double eqArea) { _eqArea = eqArea; }

    double getKElastic()const { return _kElastic; }
    void setKElastic(double kElastic) { _kElastic = kElastic; }
    
};





#endif
