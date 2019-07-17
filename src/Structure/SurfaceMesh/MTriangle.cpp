#include "Structure/SurfaceMesh/MTriangle.hpp"

#include "SysParams.h"

MTriangle::MTriangle(short membraneType) {
    if(!SysParams::Mechanics().MemCylinderVolumeK.empty())
        _kExVol = SysParams::Mechanics().MemCylinderVolumeK[membraneType];
    
}
