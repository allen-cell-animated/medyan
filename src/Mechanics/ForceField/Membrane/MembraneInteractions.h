#ifndef MEDYAN_MembraneInteractions_h
#define MEDYAN_MembraneInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Membrane;

/// Represents an internal Filament interaction
class MembraneInteractions {
    
// TODO: friend class FilamentFF;
    
protected:
    // TODO: what is it?
    /*
    /// The filament in the case of an error
    Filament* _filamentCulprit = nullptr;
    */

public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(double d) = 0; // TODO: what id d?
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};


#endif
