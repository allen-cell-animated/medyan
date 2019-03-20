//
// Created by aravind on 9/19/17.
//

#include "CUDAcommon.h"
SERLvars CUDAcommon::serlvars;
Callbacktime CUDAcommon::ctime;
Callbackcount CUDAcommon::ccount;
PolyPlusEndTemplatetime CUDAcommon::ppendtime;
#if defined(CUDAACCL) || defined(CUDATIMETRACK)
CUDAvars  CUDAcommon::cudavars;
CylCylNLvars CUDAcommon::cylcylnlvars;
SERLtime CUDAcommon::serltime;
CUDAtime CUDAcommon::cudatime;
#endif