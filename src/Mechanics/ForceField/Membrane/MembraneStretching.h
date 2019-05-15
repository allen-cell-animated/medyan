#ifndef MEDYAN_MembraneStretching_h
#define MEDYAN_MembraneStretching_h

#include "Mechanics/ForceField/Membrane/MembraneInteractions.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingHarmonic.hpp"

enum class MembraneStretchingAccumulationType {
    ByTriangle, ByVertex
};

template< MembraneStretchingAccumulationType >
class MembraneStretching: public MembraneInteractions {
private:
    MembraneStretchingHarmonic _msh;

public:
    virtual double computeEnergy(const double* coord, bool stretched) override;
    virtual void computeForces(const double* coord, double* force) override;

    virtual const string getName() { return "Membrane Stretching"; }
    
};

// Specialization declarations
template<> double MembraneStretching< MembraneStretchingAccumulationType::ByTriangle >::computeEnergy(const double* coord, bool stretched);
template<> void MembraneStretching< MembraneStretchingAccumulationType::ByTriangle >::computeForces(const double* coord, double* force);
template<> double MembraneStretching< MembraneStretchingAccumulationType::ByVertex >::computeEnergy(const double* coord, bool stretched);
template<> void MembraneStretching< MembraneStretchingAccumulationType::ByVertex >::computeForces(const double* coord, double* force);


#endif
