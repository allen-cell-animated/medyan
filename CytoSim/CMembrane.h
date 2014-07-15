//
//  CMembrane.h
//  CytoSim
//
//  Created by James Komianos on 7/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CMembrane__
#define __CytoSim__CMembrane__

#include "CFilament.h"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <iostream>


namespace chem {
    
    ///function to integrate for membrane calculation
    double gaussian(double x, void* params);
    
    /// CMembrane class updates polymerization rates based on the Brownian-Ratchet model

    /*! The CMembrane class will update all polymerization reaction rates based on
     *  the given configuration of filament ends and a membrane. The calculation is:
     *
     *
     *  1) The effective polymerization rate of a filament is equal to the "bare" rate,
     *     denoted as k0, times the probability of the gap between the filament and the
     *     membrane opening. The relation between the loading force and the rate is:
     *
     *                      k_n + = k+ * exp(-f_n * monomer size / k * T)
     *
     *     where f_n is the current force on the nth filament.
     *
     *  2) To compute the forces f_n, we must calculate the probability of finding the
     *     nth filament below the membrane (which we denote as a height of h = max_n(h_n))
     *
     *             p_n ~ integral from h-h_n to infinity ( exp(-z^2 / sigma^2) dz )
     *
     *     where sigma is the average membrane fluctuation amplitude, which we assume
     *     to be constant.
     *
     *  3) Then, the forces are computed as f_n = (p_n / p) * f, where the total force f
     *     is constant, and p is the normalization factor.
     *
     */

    class CMembrane {
        
    private:
        ///Filament and reaciton containers
        std::unordered_multimap<CFilament*, ReactionBase*> _poly_reactions; ///<map of reactions to update
        int _h_membrane = 0; /// <height of the membrane
    
        ///Membrane parameters
        double _force = 10; ///< total force on membrane (in pN)
        double _sigma = 10; ///< fluctuation amplitude (in nm)
        
        ///Integrator
        gsl_integration_workspace* _w = gsl_integration_workspace_alloc (1000);
        gsl_function _F;

    public:
        
        ///Constructor and destructor
        CMembrane() {};
        
        ~CMembrane() { gsl_integration_workspace_free(_w);}
        
        ///Add a polymerization reaction to map
        virtual void addReaction(CFilament* f, ReactionBase* r);
        
        ///Update a filament's associated polymerization reactions
        ///@note deallocates the input rxns vector
        virtual void updateFilamentReactions(CFilament* f, std::vector<ReactionBase*>* rxns);

        ///Updates the height of the membrane based on the filaments.
        ///@note see class documentation for definition of membrane height
        virtual void updateHeight();
        
        ///Updates polymerization reaction rates
        ///@note see class documentation for calculation details
        virtual void updateRates();
        
    };

}; //end namespace chem


#endif /* defined(__CytoSim__CMembrane__) */
