
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Bead_h
#define MEDYAN_Bead_h

#include <algorithm> // find
#include <vector>
#include <list>

#include "common.h"
#include "MathFunctions.h" // vec2Vector
#include "Structure/Database.h"
#include "Component.h"
#include "Composite.h"
#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "SysParams.h"
#include "Util/Math/Vec.hpp"

//FORWARD DECLARATIONS
class Compartment;
class Filament;

struct BeadData {
    using vec_type = mathfunc::Vec< 3, floatingpoint >;
    using vec_array_type = mathfunc::VecArray< 3, floatingpoint >;

    vec_array_type coords;
    vec_array_type coordsStr; // stretched coordinate
    vec_array_type forces; // currently the search dir in cg method
    vec_array_type forcesAux; // real force
    vec_array_type forcesAuxP; // prev real force

    void push_back(
        const vec_type& coord,
        const vec_type& coordStr,
        const vec_type& force,
        const vec_type& forceAux,
        const vec_type& forceAuxP
    ) {
        coords.push_back(coord);
        coordsStr.push_back(coordStr);
        forces.push_back(force);
        forcesAux.push_back(forceAux);
        forcesAuxP.push_back(forceAuxP);
    }

    void set_content(
        std::size_t pos,
        const vec_type& coord,
        const vec_type& coordStr,
        const vec_type& force,
        const vec_type& forceAux,
        const vec_type& forceAuxP
    ) {
        coords    [pos] = coord;
        coordsStr [pos] = coordStr;
        forces    [pos] = force;
        forcesAux [pos] = forceAux;
        forcesAuxP[pos] = forceAuxP;
    }

    void move_content(std::size_t from, std::size_t to) {
        coords    [to] = coords    [from];
        coordsStr [to] = coordsStr [from];
        forces    [to] = forces    [from];
        forcesAux [to] = forcesAux [from];
        forcesAuxP[to] = forcesAuxP[from];
    }

    void resize(size_t size) {
        coords    .resize(size);
        coordsStr .resize(size);
        forces    .resize(size);
        forcesAux .resize(size);
        forcesAuxP.resize(size);
    }

};

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

class Bead : public Component, public Trackable, public Movable,
    public Database< Bead, true, BeadData > {
    
public:
    using coordinate_type      = BeadData::vec_type;
    using coordinate_ref_type  = BeadData::vec_array_type::reference;
    using coordinate_cref_type = BeadData::vec_array_type::const_reference;
    using DatabaseType         = Database< Bead, true, BeadData >;

    ///@note - all vectors are in x,y,z coordinates.
    vector<floatingpoint> coordinateP; ///< Prev coordinates of bead in CG minimization

                          ///< Forces should always correspond to current coordinates.
    
    vector<floatingpoint> brforce; //boundary repulsion force
    vector<floatingpoint> pinforce;

    vector<floatingpoint> loadForcesP;
    vector<floatingpoint> loadForcesM;
    ///< The force on this bead due to an external load
    ///< This is not a vector (x,y,z) value, but a list of
    ///< force magnitudes in the direction of polymerization with
    ///< monomer increments (future values).
    ///< These are then used to propagate load forces in between
    ///< mechanical force calculations.
    ///< (Edited 20180216) Different angle between the cylinder and the
    ///< boundary would result in different effective monomer size in the
    ///< calculation of the Brownian Ratchet model. To simply computation, we
    ///< include that factor in our loadForces here. As a clarification, the
    ///< actual physical load force should not have that factor.
    
    short lfip = 0;
    short lfim = 0;  ///< Index which saves which load force to use
    
    /// The bead can be pinned to a certain position in the simulation volume.
    /// These parameters describe the pinning. Adding the Bead to the list of pinned
    /// Beads is done by a corresponding special protocol. (see executeSpecialProtocols() in Controller)
    vector<floatingpoint> pinnedPosition;
    
    bool isStatic = false;
    
    ///Main constructor
    Bead (vector<floatingpoint> v, Composite* parent, int position);
    
    ///Default constructor
    Bead(Composite* parent, int position);

    auto coordinate()    { return getDbData().coords    [getStableIndex()]; }
    auto coordinateStr() { return getDbData().coordsStr [getStableIndex()]; }
    auto force()         { return getDbData().forces    [getStableIndex()]; }
    auto forceAux()      { return getDbData().forcesAux [getStableIndex()]; }
    auto forceAuxP()     { return getDbData().forcesAuxP[getStableIndex()]; }

    auto coordinate()    const { return getDbDataConst().coords    [getStableIndex()]; }
    auto coordinateStr() const { return getDbDataConst().coordsStr [getStableIndex()]; }
    auto force()         const { return getDbDataConst().forces    [getStableIndex()]; }
    auto forceAux()      const { return getDbDataConst().forcesAux [getStableIndex()]; }
    auto forceAuxP()     const { return getDbDataConst().forcesAuxP[getStableIndex()]; }

    // Temporary compromise
    auto vcoordinate()    const { return mathfunc::vec2Vector(coordinate()   ); }
    auto vforce()         const { return mathfunc::vec2Vector(force()        ); }
    auto vforceAux()      const { return mathfunc::vec2Vector(forceAux()     ); }
    auto vforceAuxP()     const { return mathfunc::vec2Vector(forceAuxP()    ); }
    
    /// Auxilliary coordinate getter
    template< bool stretched > auto getCoordinate() const {
        return stretched ? coordinateStr() : coordinate();
    }

    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get position
    int getPosition() {return _position;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // only takes care of pinned bead removal
    virtual void addToSubSystem() override {}
    virtual void removeFromSubSystem() {
        //remove if pinned
        if(_isPinned) removeAsPinned();
    }
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Bead*>& getBeads() {
        return getElements();
    }
    
    /// Add this bead as a pinned bead
    void addAsPinned() {
        _isPinned = true;
        _pinnedBeads.push_back(this);
    }
    
    /// Remove this bead as pinned. Will remove from pinnedBeads DB
    /// @note - only usually called upon the destruction of a Bead.
    void removeAsPinned() {
        
        _isPinned = false;
        auto it = std::find(_pinnedBeads.begin(), _pinnedBeads.end(), this);
        if(it != _pinnedBeads.end()) _pinnedBeads.erase(it);
    }
    
    const vector<floatingpoint>& getPinPosition() { return pinnedPosition;}
    // Remove all pinned beads.
    void resetAllPinned() {

        _isPinned = false;
        _pinnedBeads.clear();
    }
    /// Get all pinned beads from subsystem
    static const vector<Bead*>& getPinnedBeads() {
        
        return _pinnedBeads;
    }
    
    bool isPinned() {return _isPinned;}
    
    /// Get the number of beads in this system
    static int numBeads() {
        return getElements().size();
    }
    
    /// Update the position, inherited from Movable
    virtual void updatePosition();
    
    virtual void printSelf()const;
    
    //GetType implementation just returns type of parent
    virtual int getType() {return getParent()->getType();}
    //Aravind return static
    bool getstaticstate() {return isStatic;}
    //Aravind set static
    void setstaticstate(bool index) {isStatic = index;}
    //@{
    /// Auxiliary method for CG minimization
    inline double FDotF() {
        return magnitude2(force());
    }
//    inline double FDotF() {
//        return force1[0]*force1[0] +
//        force1[1]*force1[1] +
//        force1[2]*force1[2];
//    }
    inline double FDotFA() {
        return dot(force(), forceAux());
    }
    inline double FADotFA() {
        return dot(forceAux(), forceAux());
    }
    
    inline double FADotFAP() {
        return dot(forceAux(), forceAuxP());
    }
    //Qin add brFDotbrF
    inline floatingpoint brFDotbrF() {
        return brforce[0]*brforce[0] +
        brforce[1]*brforce[1] +
        brforce[2]*brforce[2];
    }
    //add pinFDotpinF
    inline floatingpoint pinFDotpinF() {
        return pinforce[0]*pinforce[0] +
        pinforce[1]*pinforce[1] +
        pinforce[2]*pinforce[2];
    }
    //@}
    
    ///Helper functions for load forces
    
    floatingpoint getLoadForcesP();
    
    void printLoadForcesP() {
        
        cout << "loadP =";
        
        for (int i = 0; i < loadForcesP.size(); i++) {
            
            cout << " " << loadForcesP[i] << " ";
            
        }
        cout << endl;
    }
    
    floatingpoint getLoadForcesM();
 
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
    bool _isPinned = false;
    
    static std::vector<Bead*> _pinnedBeads; ///< Collection of pinned beads in SubSystem
                                         ///< (attached to some element in SubSystem)
};

#endif
