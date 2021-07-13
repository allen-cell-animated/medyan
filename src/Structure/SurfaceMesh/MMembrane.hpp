#ifndef MEDYAN_Structure_SurfaceMesh_MMembrane_hpp
#define MEDYAN_Structure_SurfaceMesh_MMembrane_hpp

#include "SysParams.h"

/******************************************************************************
Stores mechanical properties of the membrane.
******************************************************************************/

struct MMembrane {

    double kVolume  = 0; // The elastic constant of volume
    double eqVolume = 1; // The equilibrium volume

    double eqArea   = 1; // The equilibrium area, used in quadratic global stretching energy
    double kArea    = 0; // Area elasticity, used in quadratic global stretching energy

    double tension  = 0; // Membrane tension, used in linear stretching energy.
                         // Note that the tension can also be defined in quadratic energy case,
                         // but is not stored in this variable.

};



#endif
