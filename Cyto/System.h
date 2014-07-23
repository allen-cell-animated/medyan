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

#include "FilopodiaCSystem.h"
#include "ChemNRMImpl.h"
#include "CMembrane.h"
#include "MSystem.h"
#include "MNetwork.h"
using namespace chem;

template<size_t NDIM>
class FullSystem {
    
private:
    System* _mSystem;
    CSystem<NDIM>* _cSystem;
    
public:
    
    ///Constructor and destructor
    FullSystem() {}
    ~FullSystem() {}
    
    virtual void initialize(std::vector<std::vector<std::vector<double> > >& v)
    {
        ///Create chemsim, membrane
        ChemNRMImpl chem_sim_impl;
        ChemSim chem(&chem_sim_impl);
        
        CMembrane mem;
        
        ///Init filament initializer
        SimpleInitializer<1> initializer{chem, mem};

        _cSystem = new FilopodiaCSystem<NDIM>(&initializer);
        _mSystem = new System();
        _mSystem->AddNewFilaments(v,0);
        
        auto mFilamentList = _mSystem->getNetwork()[0]->getFilamentDB();
        
        int i = 0;
        for(auto &mFilament : *mFilamentList) {
            
            int length = int(mathfunc::TwoPointDistance(v[i][0], v[i][1]));
            
            auto cFilament = _cSystem->initializeCFilament(length);

            cFilament->setMFilament(mFilament);
            mFilament->setCFilament(cFilament);
            i++;
            
        }
    }
    
    virtual void print() {
        _cSystem->printFilaments();
        auto CylList = _mSystem->getCDB();
        
        for (auto it: *CylList){
            std::cout<< (*it).GetFirstBead()->coordinate[0]<< "  "
                <<(*it).GetFirstBead()->coordinate[1]<< "  "
                <<(*it).GetFirstBead()->coordinate[2]<< std::endl;}
    }
    
    
};


#endif /* defined(__CytoSim__System__) */
