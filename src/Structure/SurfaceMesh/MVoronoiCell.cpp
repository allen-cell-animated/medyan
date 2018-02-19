#include "Vertex.h"
#include "MVoronoiCell.h"

MVoronoiCell::MVoronoiCell(short membraneType) {
    
    if(!SysParams::Mechanics().MemElasticK.empty())
        _kElastic = SysParams::Mechanics().MemElasticK[membraneType];
    
    if(!SysParams::Mechanics().MemBendingK.empty())
        _kBending = SysParams::Mechanics().MemBendingK[membraneType];
    if(!SysParams::Mechanics().MemEqCurv.empty())
        _eqCurv = SysParams::Mechanics().MemEqCurv[membraneType];
}
