//
//  ConjugateGradient.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_ConjugateGradient_h
#define Cyto_ConjugateGradient_h

#include "common.h"

#include "CGFletcherRievesMethod.h"
#include "CGPolakRibiereMethod.h"

//FORWARD DECLARATIONS
class ForceFieldManager;

template <class CGType>
class ConjugateGradient : public Minimizer {
    
private:
    CGType _CGType;
    
public:
    void equlibrate(ForceFieldManager &FFM) {_CGType.minimize(FFM);}
};

#endif /* defined(__Cyto__ConjugateGradient__) */
