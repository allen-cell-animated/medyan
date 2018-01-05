#include "Triangle.h"
#include "MTriangle.h"

#include "MathFunctions.h"

using namespace mathfunc;

MTriangle::MTriangle(short membraneType) {
    if(!SysParams::Mechanics().MemElasticK.empty())
        _kElastic = SysParams::Mechanics().MemElasticK[membraneType];
    
    if(!SysParams::Mechanics().MemCylinderVolumeK.empty())
        _kExVol = SysParams::Mechanics().MemCylinderVolumeK[membraneType];
    
}
