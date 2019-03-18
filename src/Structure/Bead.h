
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

#include <vector>
#include <list>
#include "CUDAcommon.h"

#include "common.h"
#include "MathFunctions.h" // vec2Vector
#include "Structure/Database.h"
#include "Component.h"
#include "Composite.h"
#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "SysParams.h"
#include "util/math/vec.hpp"

//FORWARD DECLARATIONS
class Compartment;
class Filament;

struct BeadData {
    using vec_type = mathfunc::Vec3;
    using vec_array_type = mathfunc::VecArray< 3, double >;

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

    void pop_back() {
        coords.pop_back();
        coordsStr.pop_back();
        forces.pop_back();
        forcesAux.pop_back();
        forcesAuxP.pop_back();
    }

    void copy_from_back(std::size_t dbIndex) {
        coords[dbIndex] = coords.back();
        coordsStr[dbIndex] = coordsStr.back();
        forces[dbIndex] = forces.back();
        forcesAux[dbIndex] = forcesAux.back();
        forcesAuxP[dbIndex] = forcesAuxP.back();
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
    public Database< Bead, BeadData > {
    
public:
    using db_type = Database< Bead, BeadData >;

    ///@note - all vectors are in x,y,z coordinates.
    static bool triggercylindervectorization;
    vector<double> coordinateP; ///< Prev coordinates of bead in CG minimization

                          ///< Forces should always correspond to current coordinates.
    vector<double> force1;
    
    vector<double> brforce; //Qin boundary repulsion force
    vector<double> pinforce;

    vector<double> loadForcesP;
    vector<double> loadForcesM;
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
    vector<double> pinnedPosition;
    
    bool isStatic = false;
    
    ///Main constructor
    Bead (vector<double> v, Composite* parent, int position);
    
    ///Default constructor
    Bead(Composite* parent, int position);

    auto coordinate()    { return _dbData.coords    [_dbIndex]; }
    auto coordinateStr() { return _dbData.coordsStr [_dbIndex]; }
    auto force()         { return _dbData.forces    [_dbIndex]; }
    auto forceAux()      { return _dbData.forcesAux [_dbIndex]; }
    auto forceAuxP()     { return _dbData.forcesAuxP[_dbIndex]; }

    // Temporary compromise
    auto vcoordinate()    { return mathfunc::vec2Vector(coordinate()   ); }
    auto vcoordinateStr() { return mathfunc::vec2Vector(coordinateStr()); }
    auto vforce()         { return mathfunc::vec2Vector(force()        ); }
    auto vforceAux()      { return mathfunc::vec2Vector(forceAux()     ); }
    auto vforceAuxP()     { return mathfunc::vec2Vector(forceAuxP()    ); }
    
    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get position
    int getPosition() {return _position;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() override {}
    virtual void removeFromSubSystem() {
        /* Haoran 03/17/2019
//        std::cout<<"removing bead with bindex "<<_dbIndex<<endl;
        //Reset in bead coordinate vector and add _dbIndex to the list of removedbindex.
        removedbindex.push_back(_dbIndex);
        resetcoordinates();
        */
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
        _pinnedBeads.addElement(this);
    }
    
    /// Remove this bead as pinned. Will remove from pinnedBeads DB
    /// @note - only usually called upon the destruction of a Bead.
    void removeAsPinned() {
        
        _isPinned = false;
        _pinnedBeads.removeElement(this);
    }
    
    const vector<double>& getPinPosition() { return pinnedPosition;}

    //Qin
    // Remove all pinned beads.
    void resetAllPinned() {

        _isPinned = false;
        _pinnedBeads.clearElements();
    }
    /// Get all pinned beads from subsystem
    static const vector<Bead*>& getPinnedBeads() {
        
        return _pinnedBeads.getElements();
    }
    
    bool isPinned() {return _isPinned;}
    
    /// Get the number of beads in this system
    static int numBeads() {
        return getElements().size();
    }
    
    /// Update the position, inherited from Movable
    virtual void updatePosition();
    
    virtual void printSelf();
    
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
    inline double brFDotbrF() {
        return brforce[0]*brforce[0] +
        brforce[1]*brforce[1] +
        brforce[2]*brforce[2];
    }
    //Qin add pinFDotpinF
    inline double pinFDotpinF() {
        return pinforce[0]*pinforce[0] +
        pinforce[1]*pinforce[1] +
        pinforce[2]*pinforce[2];
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
    
    bool _isPinned = false;
    
    static OldDatabase<Bead*> _pinnedBeads; ///< Collection of pinned beads in SubSystem
                                         ///< (attached to some element in SubSystem)
    //Vectorize beads so the coordinates are all available in a single array.
    //@{
    static int vectormaxsize;//maximum number of beads that can be appended without
    // revectorization
    [[deprecated]] static int Nbeads;//Total number of beads in the system
    static vector<int> removedbindex;//stores the bead indices that have been freed
    // through depolymerization/ destruction reactions.
    /* Haoran 03/17/2019
    static void revectorizeifneeded(){
        int newsize = vectormaxsize;
        //if the maximum bead index is very close to the vector size
        if(vectormaxsize - maxbindex <= bead_cache/10 )
            //new size will be increased by bead_cache
            newsize = (int(Nbeads/bead_cache)+2)*bead_cache;
        //if we have removed bead_cache number of beads from the system
        if(removedbindex.size() >= bead_cache)
            //we can revectorize with a smaller size.
            newsize = (int(Nbeads/bead_cache)+1)*bead_cache;
        //set parameters and revectorize
        if(newsize != vectormaxsize){
            double *coord = CUDAcommon::serlvars.coord;
            delete[] coord;
            double *newcoord = new double[3 * newsize];
            CUDAcommon::serlvars.coord = newcoord;
            revectorize(newcoord);
            //copyvector(newcoord, coord);
            vectormaxsize = newsize;
            //cylinder structure needs to be revecotrized as well.
            triggercylindervectorization = true;
        }
    }
    static void revectorize(double* coord){
        //set contiguous bindices and set coordinates.
        int idx = 0;
        for(auto b:_beads.getElements()){
            int index = 3 * idx;
            coord[index] = b->coordinate[0];
            coord[index + 1] = b->coordinate[1];
            coord[index + 2] = b->coordinate[2];
            b->getDbIndex() = idx;
            idx++;
        }
        Nbeads =_beads.getElements().size();
        maxbindex = _beads.getElements().size();
        removedbindex.clear();
    }*/
    /* Haoran 03/17/2019
    //copy coodinates of this bead to the appropriate spot in coord vector.
    void  copycoordinatestovector() {
//        if(!triggercylindervectorization) {
            CUDAcommon::serlvars.coord[3 * _dbIndex] = coordinate[0];
            CUDAcommon::serlvars.coord[3 * _dbIndex + 1] = coordinate[1];
            CUDAcommon::serlvars.coord[3 * _dbIndex + 2] = coordinate[2];
//        }
    }*/
    void resetcoordinates() {
        CUDAcommon::serlvars.coord[3 * _dbIndex] = -1.0;
        CUDAcommon::serlvars.coord[3 * _dbIndex + 1] = -1.0;
        CUDAcommon::serlvars.coord[3 * _dbIndex + 2] = -1.0;
    }
    //@}
};


#endif
