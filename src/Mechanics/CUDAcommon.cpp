//
// Created by aravind on 9/19/17.
//

#include "CUDAcommon.h"
Callbacktime CUDAcommon::ctime;
Callbackcount CUDAcommon::ccount;
PolyPlusEndTemplatetime CUDAcommon::ppendtime;
timeminimization CUDAcommon::tmin;
motorwalking CUDAcommon::mwalk;
chemdetails CUDAcommon::cdetails;
#if defined(CUDAACCL) || defined(CUDATIMETRACK)
CUDAvars  CUDAcommon::cudavars;
CylCylNLvars CUDAcommon::cylcylnlvars;
SERLtime CUDAcommon::serltime;
CUDAtime CUDAcommon::cudatime;
#endif