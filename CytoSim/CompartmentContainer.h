//
//  CompartmentContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 9/5/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CompartmentContainer__
#define __CytoSim__CompartmentContainer__

#include <iostream>
#include "Compartment.h"

///Key to access instance of CompartmentGrid
class CompartmentGridKey {friend class SubSystem;
                          friend class SimpleInitializerImpl;
                          friend class GController;
                          CompartmentGridKey(){} ~CompartmentGridKey(){} };


/// CompartmentGrid class is a simple n-dimensional grid of CompartmentSpatial objects (singleton)

/*!
 *  The CompartmentGrid class is a singleton grid of CompartmentSpatial objects, each of them seperately
 *  holding internal and diffusion reactions, species information, as well as spatial information.
 *  This class is n-dimensional, and the dimension is specified at runtime.
 *
 *  All compartments within CompartmentGrid are indexed, and this class is responsible for
 *  assignment of compartment neighbors upon initialization.
 *
 *  The _prototype_compartment is used such that all reactions and species can be added to 
 *  the prototype, and this configuration can be copied to all compartments in the grid.
 *  An example of such is below:
 *
 *  @code
 *  CompartmentGrid::setInstance(CompartmentGridKey(), {50,50,50});
 *  Compartment &Cproto = CompartmentGrid::Instance()->getProtoCompartment();
 *  Species *M1 = Cproto.addSpecies("Myosin",1U);
 *  Cproto.setDiffusionRate(M1,2000);
 *  Species *M2 = Cproto.addSpecies("Fascin",6U);
 *  Cproto.setDiffusionRate(M2,2000);
 *
 *  vector<float> sides{100.0,100.0,100.0};
 *  Cproto.setSides(sides.begin());
 *
 *  Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
 *  Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
 *  CompartmentGrid::Instance()->initialize();
 *  @endcode
 *
 *  @note the protocompartment's side length must be specified so that all compartments will
 *  be initialized to the same side length. The protocomparment must be set up BEFORE initalize()
 *  is called.
 */
class CompartmentGrid : public Composite {
private:
    Compartment _prototype_compartment; ///< prototype compartment, to be configured before initialization
    bool _is_initialized; /// grid is initalized or not
    int _numCompartments;
    
    static CompartmentGrid* _instance; ///singleton instance
    
    //private constructor
    CompartmentGrid(int numCompartments) : _numCompartments(numCompartments) {}
public:
    
    /// Copying is not allowed
    CompartmentGrid(const CompartmentGrid &rhs) = delete;
    
    /// Assignment is not allowed
    CompartmentGrid& operator=(CompartmentGrid &rhs) = delete;
    
    ///set instance of grid (should only be done at beginning of program)
    static void setInstance(CompartmentGridKey k, int numCompartments);
    
    ///Get instance of grid
    static CompartmentGrid* Instance(CompartmentGridKey k);
    
    /// Get name of this compartment grid
    virtual std::string getFullName() const {return std::string("CompartmentGrid");};
    
    /// Initialize the compartment, which copies all species and reactions of the protocompartment into the
    /// compartments in the grid. Also generates neighboring connections for all compartments as well as
    /// initializes spatial coordinates of all compartments.
    /// @note - _prototype_compartment must be configured before this is called
    void initialize(int length)
    {
        if(_is_initialized)
            throw std::runtime_error("CompartmentGrid::initialize() should be called only once");
        
        for(size_t i=0; i<length; ++i)
        {
            addChild(std::unique_ptr<Component>(new Compartment()));
        }
        
        _is_initialized=true;
    }
    
    ///Activate all compartments
    void activateAll()
    {
        for (auto &c : children())
            static_cast<Compartment*>(c.get())->activate();
        
    }
    
    /// Get the protocompartment from this grid, in order to configure and then initialize
    Compartment& getProtoCompartment() {return _prototype_compartment;}
    const Compartment& getProtoCompartment() const {return _prototype_compartment;}
    
    
    /// Add reactions to all compartments in the grid
    /// @param - chem, a ChemSim object that controls the reaction algorithm
    virtual void addChemSimReactions()
    {
        for(auto &c : children())
        {
            Compartment*C = static_cast<Compartment*>(c.get());
            C->addChemSimReactions();
        }
    }
    
    /// Print properties of this grid
    virtual void printSelf()
    {
        std::cout << getFullName() << std::endl;
        std::cout << "Number of Compartment objects: " << numberOfChildren() << std::endl;
        for(auto &c : children())
            c->printSelf();
    }

};


#endif /* defined(__CytoSim__CompartmentContainer__) */
