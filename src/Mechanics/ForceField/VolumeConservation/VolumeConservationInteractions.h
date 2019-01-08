#ifndef MEDYAN_VolumeConservationInteractions_h
#define MEDYAN_VolumeConservationInteractions_h

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
    virtual double computeEnergy(bool stretched) = 0; // d is the stretching parameter along the force
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the name of this interaction
    virtual string getName()const = 0;
};


#endif
