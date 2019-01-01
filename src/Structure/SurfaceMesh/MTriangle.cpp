#include "Triangle.h"
#include "MTriangle.h"

MTriangle::MTriangle(short membraneType) {
    if(!SysParams::Mechanics().MemCylinderVolumeK.empty())
        _kExVol = SysParams::Mechanics().MemCylinderVolumeK[membraneType];
    
}
