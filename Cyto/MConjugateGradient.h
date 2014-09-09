//
//  MConjugateGradient.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MConjugateGradient_h
#define Cyto_MConjugateGradient_h

class MController;
class Minimize;

template <class CGType>
class ConjugateGradient : private Minimize
{
    
private:
    CGType _CGType;
    
    
public:
    
    void Equlibrate(MController*){_CGType.Minimize(MController*);}
    
};

#endif
