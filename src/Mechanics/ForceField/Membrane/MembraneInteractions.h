#ifndef MEDYAN_MembraneInteractions_h
#define MEDYAN_MembraneInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Membrane;

/// Represents an internal Filament interaction
class MembraneInteractions {
    
friend class MembraneFF;
    
protected:
    /// The membrane in the case of an error
    Membrane* _membraneCulprit = nullptr;

public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(const double* coord, bool stretched) = 0; // d is the stretching parameter along the force
    /// Compute forces of this interaction
    virtual void computeForces(const double* coord, double* force) = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};


#endif
