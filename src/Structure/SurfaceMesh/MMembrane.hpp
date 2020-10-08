#ifndef MEDYAN_Structure_SurfaceMesh_MMembrane_hpp
#define MEDYAN_Structure_SurfaceMesh_MMembrane_hpp

#include "SysParams.h"

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

struct MMembrane {

    double eqVolume = 1; // The equilibrium volume

    double eqArea   = 1; // The equilibrium area, used in quadratic stretching energy
    double kElastic = 0; // Area elasticity, used in quadratic stretching energy

    double tension  = 0; // Membrane tension, used in linear stretching energy.
                         // Note that the tension can also be defined in quadratic energy case,
                         // but is not stored in this variable.

    MMembrane(short membraneType) {
        if(!SysParams::Mechanics().memAreaK.empty())
            kElastic = SysParams::Mechanics().memAreaK[membraneType];
        if(!SysParams::Mechanics().MemTension.empty())
            tension = SysParams::Mechanics().MemTension[membraneType];
    }

    double getEqVolume()const { return eqVolume; }

    double getEqArea()const { return eqArea; }

    double getKElastic()const { return kElastic; }

};





#endif
