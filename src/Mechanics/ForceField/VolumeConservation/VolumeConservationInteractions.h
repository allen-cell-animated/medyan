#ifndef MEDYAN_VolumeConservationInteractions_h
#define MEDYAN_VolumeConservationInteractions_h

#include <string>

//FORWARD DECLARATIONS
class Membrane;

/// Represents an internal Filament interaction
class VolumeConservationInteractions {
    
    friend class VolumeConservationFF;
    
protected:
    /// The membrane in the case of an error
    Membrane* _membraneCulprit = nullptr;

public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(const double* coord, bool stretched) = 0; // d is the stretching parameter along the force
    /// Compute forces of this interaction
    virtual void computeForces(const double* coord, double* force) = 0;
    
    /// Get the name of this interaction
    virtual std::string getName()const = 0;
};


#endif
