#include "Structure/SurfaceMesh/MTriangle.hpp"

#include "SysParams.h"

MTriangle::MTriangle(short membraneType) {
    if(!SysParams::Mechanics().MemBeadVolumeK.empty())
        _kExVol = SysParams::Mechanics().MemBeadVolumeK[membraneType];
    
}
