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
        //int maxCompSize;
        size_t NGRID = 5;
        CompartmentGrid<NDIM>* grid;
        
        switch (NDIM) {
            case 1:
                grid = new CompartmentGrid<NDIM>({NGRID});
                break;
            case 2:
                grid = new CompartmentGrid<NDIM>({NGRID,NGRID});
                break;
            default:
                grid = new CompartmentGrid<NDIM>({NGRID,NGRID,NGRID});
        }
        
        
        CompartmentSpatial<NDIM> &Cproto = grid->getProtoCompartment();
        
        Species *M1 = Cproto.addSpecies("Actin",10U);
        Species *M2 = Cproto.addSpecies("Capping",0U);
        Species *M3 = Cproto.addSpecies("X-Formin",5U);
        Species *M4 = Cproto.addSpecies("Myosin", 5U);
        Cproto.setDiffusionRate(M1,2000);
        Cproto.setDiffusionRate(M2,2000);
        Cproto.setDiffusionRate(M3,2000);
        Cproto.setDiffusionRate(M4,2000);

        ///Set side length
        std::vector<float> sides{100.0};
        Cproto.setSides(sides.begin());
        
        ///init grid
        grid->initialize();
        
        ///Create chemsim and init
        ChemNRMImpl chem_sim_impl;
        ChemSim chem(&chem_sim_impl);
        grid->addChemSimReactions(chem);
        
        CMembrane mem;
        
        ///Init filament initializer
        SimpleInitializer<1> initializer{chem, mem};

        _cSystem = new FilopodiaCSystem<NDIM>(grid,&initializer);
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
