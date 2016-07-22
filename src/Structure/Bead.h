
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Bead_h
#define MEDYAN_Bead_h

#include <vector>
#include <list>

#include "common.h"

#include "Database.h"
#include "Composite.h"
#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "SysParams.h"

//FORWARD DECLARATIONS
class Compartment;
class Filament;

/// Represents a single coordinate between [Cylinders](@ref Cylinder), and holds forces
/// needed for mechanical equilibration.
/*!
 *  Beads are the "hinges" between [Cylinders](@ref Cylinder). In the minimization 
 *  algorithms, beads are moved corresponding to external forces, for example, Filament 
 *  stretching and bending. The bead class contains currernt coordinates and forces, and 
 *  has functions to calculate dot products for the minimization algorithms.
 *
 *  Extending the Movable class, the positions of all instances can 
 *  be updated by the SubSystem.
 *
 *  Extending the DynamicNeighbor class, all instances can be kept in 
 *  [NeighborLists](@ref NeighborList).
 */

class Bead : public Component, public Trackable, public Movable{
    
public:
    ///@note - all vectors are in x,y,z coordinates.
    
    vector<double> coordinate;  ///< Coordinates of the bead
    vector<double> coordinateP; ///< Prev coordinates of bead in CG minimization
    vector<double> coordinateB; ///< Prev coordinate of bead before CG minimization
    
	vector<double> force; ///< Forces based on curent coordinates.
                          ///< Forces should always correspond to current coordinates.
    vector<double> forceAux;  ///< An auxiliary field needed during CG minimization.
    vector<double> forceAuxP; ///< An auxiliary field needed during CG minimization.
    
    vector<double> loadForcesP;
    vector<double> loadForcesM;
    ///< The force on this bead due to an external load
    ///< This is not a vector (x,y,z) value, but a list of
    ///< force magnitudes in the direction of polymerization with
    ///< monomer increments (future values).
    ///< These are then used to propagate load forces in between
    ///< mechanical force calculations.
    
    short lfip = 0; 
    short lfim = 0;  ///< Index which saves which load force to use
    
    /// The bead can be pinned to a certain position in the simulation volume.
    /// These parameters describe the pinning. Adding the Bead to the list of pinned
    /// Beads is done by a corresponding special protocol. (see executeSpecialProtocols() in Controller)
    vector<double> pinnedPosition;
    
    ///Main constructor
    Bead (vector<double> v, Composite* parent, int position);
    
    ///Default constructor
    Bead(Composite* parent, int position);
    ~Bead() {}
    
    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get position
    int getPosition() {return _position;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _beads.addElement(this);}
    virtual void removeFromSubSystem() {_beads.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Bead*>& getBeads() {
        return _beads.getElements();
    }
    
    /// Add this bead as a pinned bead
    void addAsPinned() {
        _pinnedBeads.addElement(this);
    }
    
    /// Get all pinned beads from subsystem
    static const vector<Bead*>& getPinnedBeads() {
        
        return _pinnedBeads.getElements();
    }
    
    /// Get the number of beads in this system
    static int numBeads() {
        return _beads.countElements();
    }
    
    /// Update the position, inherited from Movable
    virtual void updatePosition();
    
    virtual void printSelf();
    
    //GetType implementation just returns type of parent
    virtual int getType() {return getParent()->getType();}
    //@{
    /// Auxiliary method for CG minimization
    inline double FDotF() {
        return force[0]*force[0] +
               force[1]*force[1] +
               force[2]*force[2];
    }
    inline double FDotFA() {
        return force[0]*forceAux[0] +
               force[1]*forceAux[1] +
               force[2]*forceAux[2];
    }
    
    inline double FADotFA() {
        return forceAux[0]*forceAux[0] +
               forceAux[1]*forceAux[1] +
               forceAux[2]*forceAux[2];
    }
    
    inline double FADotFAP() {
        return forceAux[0]*forceAuxP[0] +
               forceAux[1]*forceAuxP[1] +
               forceAux[2]*forceAuxP[2];
    }
    //@}
    
    ///Helper functions for load forces
    
    double getLoadForcesP();
    
    void printLoadForcesP() {
        
        cout << "loadP =";
        
        for (int i = 0; i < loadForcesP.size(); i++) {
            
            cout << " " << loadForcesP[i] << " ";
            
        }
        cout << endl;
    }
    
    double getLoadForcesM();
 
    void printLoadForcesM()  {
        
        cout << "loadM =";
        
        for (int i = 0; i < loadForcesM.size(); i++) {
            
            cout << " " << loadForcesM[i] << " ";
            
        }
        cout << endl;
    }
    
private:
    Compartment* _compartment = nullptr; ///< Pointer to the compartment that this bead is in
    
    int _position;     ///< Position on structure
    float _birthTime;  ///< Time of birth
    
    static Database<Bead*> _beads; ///< Collection of beads in SubSystem
    static Database<Bead*> _pinnedBeads; ///< Collection of pinned beads in SubSystem
                                         ///< (attached to some element in SubSystem)
};


#endif
