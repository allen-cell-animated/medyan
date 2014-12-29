
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Bead_h
#define M3SYM_Bead_h

#include <vector>
#include <list>

#include "common.h"

#include "BeadDB.h"

#include "Component.h"
#include "Neighbor.h"
#include "Movable.h"

//FORWARD DECLARATIONS
class Compartment;

/// Represents a single coordinate between [Cylinders](@ref Cylinder), and holds forces
/// needed for mechanical equilibration.
/*!
 *  Beads are the "hinges" between [Cylinders](@ref Cylinder). In the minimization 
 *  algorithms, beads are moved corresponding to external forces, for example, Filament 
 *  stretching and bending. The bead class contains currernt coordinates and forces, and 
 *  has functions to calculate dot products for the minimization algorithms.
 */

class Bead : public Component, public Neighbor, public Movable {
public:
    vector<double> coordinate; ///< Coordinates of the bead
    vector<double> coordinateAux; ///< An auxiliary coordinate field needed
                                  ///< during CG minimization
	vector<double> force; ///< Forces based on curent coordinates.
                          ///< Forces should always correspond to current coordinates.
    vector<double> forceAux; ///< An auxiliary field needed during CG minimization.
    
    double loadForce = 0.0; ///< The force on this bead due to an external load
                            ///< Usually a boundary element
    
    ///Main constructor
    Bead (vector<double> v, int positionFilament);
    
    ///Default constructor
    Bead(int positionFilament) :
        _positionFilament(positionFilament), coordinate (3, 0),
        coordinateAux(3, 0), force(3, 0), forceAux(3, 0) {}
    
    ~Bead();
    
    //@{
    /// Auxiliary method for CG minimization
    inline double calcForceSquare() {
        return force[0]*force[0] +
               force[1]*force[1] +
               force[2]*force[2];
    }
    
    inline double calcForceAuxSquare() {
        return forceAux[0]*forceAux[0] +
               forceAux[1]*forceAux[1] +
               forceAux[2]*forceAux[2];
    }
    
    inline double calcDotForceProduct() {
        return force[0]*forceAux[0] +
               force[1]*forceAux[1] +
               force[2]*forceAux[2];
    }
    //@}
    
    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Update the position
    virtual void updatePosition();
    
    /// Set position on the local Filament
    void setPositionFilament(int positionFilament) {
        _positionFilament = positionFilament;
    }
    /// Get position on the local Filament
    int getPositionFilament() {return _positionFilament;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}

    
private:
    Compartment* _compartment = nullptr;
        ///< Pointer to the compartment that this bead is in
    
    int _positionFilament; ///< Position on Filament
    float _birthTime; ///< Time of birth
};


#endif
