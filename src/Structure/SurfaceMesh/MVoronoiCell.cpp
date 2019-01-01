#include "Structure/SurfaceMesh/MVoronoiCell.h"

#include "SysParams.h"

MVoronoiCell::MVoronoiCell(short membraneType) {
    
    if(!SysParams::Mechanics().MemBendingK.empty())
        _kBending = SysParams::Mechanics().MemBendingK[membraneType];
    if(!SysParams::Mechanics().MemEqCurv.empty())
        _eqCurv = SysParams::Mechanics().MemEqCurv[membraneType];
}
