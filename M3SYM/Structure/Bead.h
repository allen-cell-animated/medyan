
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Bead_h
#define M3SYM_Bead_h

#include <vector>
#include <list>

#include "common.h"

#include "Database.h"
#include "Component.h"
#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"

//FORWARD DECLARATIONS
class Compartment;

/// Represents a single coordinate between [Cylinders](@ref Cylinder), and holds forces
/// needed for mechanical equilibration.
/*!
 *  Beads are the "hinges" between [Cylinders](@ref Cylinder). In the minimization 
 *  algorithms, beads are moved corresponding to external forces, for example, Filament 
 *  stretching and bending. The bead class contains currernt coordinates and forces, and 
 *  has functions to calculate dot products for the minimization algorithms.
 *  Extending the Movable and Trackable class, all instances are kept, can be moved,
 *  and easily accessed by the SubSystem.
 */

class Bead : public Component, public Trackable, public Movable, public DynamicNeighbor {
public:
    vector<double> coordinate; ///< Coordinates of the bead
    vector<double> coordinateP; ///< Prev coordinates of bead
    
	vector<double> force; ///< Forces based on curent coordinates.
                          ///< Forces should always correspond to current coordinates.
    vector<double> forceAux;  ///< An auxiliary field needed during CG minimization.
    vector<double> forceAuxP; ///< An auxiliary field needed during CG minimization.
    
    double loadForce = 0.0; ///< The force on this bead due to an external load
                            ///< Usually a boundary element
    
    ///Main constructor
    Bead (vector<double> v, int positionFilament);
    
    ///Default constructor
    Bead(int positionFilament)
        : coordinate (3, 0), force(3, 0), forceAux(3, 0),
          _positionFilament(positionFilament) {}
    
    ~Bead();
    
    //@{
    /// Auxiliary method for CG minimization
    inline double FDotF() {
        return force[0]*force[0]+
               force[1]*force[1]+
               force[2]*force[2];
    }
    inline double FDotFA() {
        return force[0]*forceAux[0]+
               force[1]*forceAux[1]+
               force[2]*forceAux[2];
    }
    
    inline double FADotFA() {
        return forceAux[0]*forceAux[0]+
               forceAux[1]*forceAux[1]+
               forceAux[2]*forceAux[2];
    }
    
    inline double FADotFAP() {
        return forceAux[0]*forceAuxP[0]+
               forceAux[1]*forceAuxP[1]+
               forceAux[2]*forceAuxP[2];
    }
    //@}
    
    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Set position on the local Filament
    void setPositionFilament(int positionFilament) {
        _positionFilament = positionFilament;
    }
    /// Get position on the local Filament
    int getPositionFilament() {return _positionFilament;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _beads.addElement(this);}
    virtual void removeFromSubSystem() {_beads.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const unordered_set<Bead*>& getBeads() {
        return _beads.getElements();
    }
    /// Get the number of beads in this system
    static int numBeads() {
        return _beads.countElements();
    }
    
    /// Update the position, inherited from Movable
    virtual void updatePosition();
    
    
private:
    Compartment* _compartment = nullptr;
        ///< Pointer to the compartment that this bead is in
    
    int _positionFilament; ///< Position on Filament
    float _birthTime; ///< Time of birth
    
    static Database<Bead*> _beads; ///< Collection of beads in SubSystem
};


#endif
