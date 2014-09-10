//
//  MConjugateGradient.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MConjugateGradient_h
#define Cyto_MConjugateGradient_h

#include "MCGFletcherRievesMethod.h"
#include "MCGPolakRibiereMethod.h"

class ForceFieldManager;
class Minimizer;


template <class CGType>
class ConjugateGradient : public Minimizer
{
    
private:
    CGType _CGType;
    
public:
    void Equlibrate(ForceFieldManager &FFM) {_CGType.Minimize(FFM);}
};

#endif
