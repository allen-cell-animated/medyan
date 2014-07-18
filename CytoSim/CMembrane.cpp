//
//  CMembrane.cpp
//  CytoSim
//
//  Created by James Komianos on 7/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CMembrane.h"

namespace chem {
    
    ///function to integrate for membrane calculation
    double gaussian(double x, void* params) {
        double sigma = *(double *) params;
        return exp(- pow(x,2) / pow(sigma,2));
    }
    
    ///Add a polymerization reaction to map
    void CMembrane::addReaction(CFilament* f, ReactionBase* r) {
        _poly_reactions.insert(std::pair<CFilament*, ReactionBase*>(f,r));
        
    }
    
    
    ///Update a filament's associated polymerization reactions
    ///@note deallocates the input rxns vector
    void CMembrane::updateFilamentReactions(CFilament* f, std::vector<ReactionBase*>* rxns)
    {
        ///Erase old rxns, add new ones
        _poly_reactions.erase(f);
        
        for(auto it = rxns->begin(); it != rxns->end(); it++)
            CMembrane::addReaction(f, (*it));
        //clean up
        delete rxns;
    }

    ///Updates the height of the membrane based on the filaments.
    ///@note see class documentation for definition of membrane height
    void CMembrane::updateHeight() {
        
        for(auto it = _poly_reactions.begin(); it != _poly_reactions.end(); it++)
        {
            int filamentHeight = (*it).first->length() * monomer_size;
            if((*it).first->length() > _h_membrane) _h_membrane = filamentHeight;
        }
    }
    
    ///Updates polymerization reaction rates
    ///@note see class documentation for calculation details
    void CMembrane::updateRates() {
        
        //Vectors for calculations
        auto fValues = std::vector<double>();
        auto pValues = std::vector<double>();
        double sumP = 0;
        
        
        //Set up integration function
        _F.function = &gaussian;
        _F.params = &_sigma;
        
        //Calculate pN, store in array
        for(auto it = _poly_reactions.begin(); it != _poly_reactions.end(); it++)
        {
            int h = (*it).first->length() * monomer_size;
            
            double p, error;
            gsl_integration_qags(&_F, h - _h_membrane, 100, 1e-6, 1e-6, 1000, _w, &p, &error);
            
            pValues.push_back(p);
            sumP += p;
        }
        
        ///Calculate forces
        for(auto it = pValues.begin(); it != pValues.end(); it++)
        {
            double pN = (*it);
            fValues.push_back( (pN / sumP) * _force);
        }
        
        ///Update rates
        int index = 0;
        for(auto it = _poly_reactions.begin(); it != _poly_reactions.end(); it++)
        {
            ReactionBase* r = (*it).second;
            double newRate = r->getBareRate() * exp(- (fValues[index++] * monomer_size) / kT );
            
            r->setRate(newRate);
            r->activateReaction();
        }
    }

}

