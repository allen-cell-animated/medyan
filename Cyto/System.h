//
//  System.h
//  CytoSim
//
//  Created by James Komianos on 7/21/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__System__
#define __CytoSim__System__

#include <iostream>

#include "CSystem.h"
#include "MSystem.h"

class System {
    
private:
    MSystem* _mSystem;
    CSystem* _cSystem;
    
public:
    
    System() {}
    
    ~System() {}
    
    virtual void initialize(std::vector<std::vector<std::vector<double> > > v, int ndim )
    {
        int maxCompSize;
        
        if(ndim == 1)
            _cSystem = new CSystem<ndim>({maxCompSize});
        else if(ndim == 2)
            _cSystem = new CSystem<ndim>({maxCompSize, maxCompSize});
        else
            _cSystem = new CSystem<ndim>({maxCompSize, maxCompSize, maxCompSize});
        
        _mSystem = new MSystem();
        
        _mSystem->AddNewFilaments(v,0);
        
        auto mFilamentList = _mSystem->getNetwork()[0]->getFilamentDB();
        
        int i = 0;
        for(auto &mFilament : mFilamentList) {
            
            int length = int(mathfunc::TwoPointDistance(v[i][0], v[i][1]));
            
            
            i++;
            
        }
    }
    
    
    
    
};


#endif /* defined(__CytoSim__System__) */
