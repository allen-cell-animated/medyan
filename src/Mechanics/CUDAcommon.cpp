//
// Created by aravind on 9/19/17.
//

#include "CUDAcommon.h"
SERLvars CUDAcommon::serlvars;
Callbacktime CUDAcommon::ctime;
Callbackcount CUDAcommon::ccount;
PolyPlusEndTemplatetime CUDAcommon::ppendtime;
timeminimization CUDAcommon::tmin;
motorwalking CUDAcommon::mwalk;
#if defined(CUDAACCL) || defined(CUDATIMETRACK)
CUDAvars  CUDAcommon::cudavars;
CylCylNLvars CUDAcommon::cylcylnlvars;
SERLtime CUDAcommon::serltime;
CUDAtime CUDAcommon::cudatime;
#endif