
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
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

#include "Database.h"
#include "Component.h"
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
    static bool triggercylindervectorization;
    vector<floatingpoint> coordinate;  ///< Coordinates of the bead
    vector<floatingpoint> coordinateP; ///< Prev coordinates of bead in CG minimization
    int _ID; ///<Bead IDs
    int _dbIndex =  -1; ///<Position in database vector

	vector<floatingpoint> force; ///< Forces based on curent coordinates.
                          ///< Forces should always correspond to current coordinates.
    vector<floatingpoint> forceAux;  ///< An auxiliary field needed during CG minimization.
    vector<floatingpoint> forceAuxP; ///< An auxiliary field needed during CG minimization.
    
    vector<floatingpoint> brforce; //Qin boundary repulsion force
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
    
    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get position
    int getPosition() {return _position;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _beads.addElement(this);}
    virtual void removeFromSubSystem() {
        //Reset in bead coordinate vector and add _dbIndex to the list of removedbindex.
        removedbindex.push_back(_dbIndex);
        resetcoordinates();
        //remove from database
        _beads.removeElement(this);
        //remove if pinned
        if(_isPinned) removeAsPinned();
        Nbeads = _beads.getElements().size();
    }
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Bead*>& getBeads() {
        return _beads.getElements();
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
    
    const vector<floatingpoint>& getPinPosition() { return pinnedPosition;}
    //Qin
    // Remove all pinned beads.
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
        return _beads.countElements();
    }
    
    /// Update the position, inherited from Movable
    virtual void updatePosition();
    
    virtual void printSelf();
    
    //GetType implementation just returns type of parent
    virtual int getType() {return getParent()->getType();}
    //Aravind get ID
    int getID() {return _ID;}
    //Aravind return static
    bool getstaticstate() {return isStatic;}
    //Aravind set static
    void setstaticstate(bool index) {isStatic = index;}
    //@{
    /// Auxiliary method for CG minimization
    inline floatingpoint FDotF() {
        return force[0]*force[0] +
               force[1]*force[1] +
               force[2]*force[2];
    }
    inline floatingpoint FDotFA() {
        return force[0]*forceAux[0] +
               force[1]*forceAux[1] +
               force[2]*forceAux[2];
    }
    
    inline floatingpoint FADotFA() {
        return forceAux[0]*forceAux[0] +
               forceAux[1]*forceAux[1] +
               forceAux[2]*forceAux[2];
    }
    
    inline floatingpoint FADotFAP() {
        return forceAux[0]*forceAuxP[0] +
               forceAux[1]*forceAuxP[1] +
               forceAux[2]*forceAuxP[2];
    }
    //Qin add brFDotbrF
    inline floatingpoint brFDotbrF() {
        return brforce[0]*brforce[0] +
        brforce[1]*brforce[1] +
        brforce[2]*brforce[2];
    }
    //Qin add pinFDotpinF
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

    static int getmaxbindex(){
        return maxbindex;
    }

	// through depolymerization/ destruction reactions.
	static void revectorizeifneeded(){
		//Run the special protocol during chemistry, the regular otherwise.
		if(SysParams::DURINGCHEMISTRY)
			appendrevectorizeifneeded();
		else {
			int newsize = vectormaxsize;
			//if the maximum bead index is very close to the vector size
			if (vectormaxsize - maxbindex <= bead_cache / 10)
				//new size will be increased by bead_cache
				newsize = (int(Nbeads / bead_cache) + 2) * bead_cache;
			//if we have removed bead_cache number of beads from the system
			if (removedbindex.size() >= bead_cache)
				//we can revectorize with a smaller size.
				newsize = (int(Nbeads / bead_cache) + 1) * bead_cache;
			//set parameters and revectorize
			if (newsize != vectormaxsize) {
//				cout<<"vectorize bead"<<endl;
				floatingpoint *coord = CUDAcommon::serlvars.coord;
				delete[] coord;
				floatingpoint *newcoord = new floatingpoint[3 * newsize];
				CUDAcommon::serlvars.coord = newcoord;
				revectorize(newcoord);
				//copyvector(newcoord, coord);
				vectormaxsize = newsize;
				//cylinder structure needs to be revecotrized as well.
				triggercylindervectorization = true;
			}
		}
/*		cout<<"Printing bead data triggercylindervectorization "
		""<<triggercylindervectorization<<endl;
		for(auto b:_beads.getElements()){
			cout<<"Bead ID "<<b->getID()<<" dbIndex "<<b->_dbIndex<<endl;
		}
		cout<<"removedbindex "<<removedbindex.size()<<endl;
		cout<<"------------------------------Bead"<<endl;*/
	}
	static void printBeaddata(){
		cout<<"Printing bead data "<<endl;
		for(auto b:_beads.getElements()){
			cout<<"Bead ID "<<b->getID()<<" dbIndex "<<b->_dbIndex<<endl;
		}
		cout<<"removedbindex "<<removedbindex.size()<<endl;
		cout<<"------------------------------Bead"<<endl;
    }
private:
    Compartment* _compartment = nullptr; ///< Pointer to the compartment that this bead is in
    
    int _position;     ///< Position on structure
    float _birthTime;  ///< Time of birth
    
    bool _isPinned = false;
    
    static Database<Bead*> _beads; ///< Collection of beads in SubSystem
    static Database<Bead*> _pinnedBeads; ///< Collection of pinned beads in SubSystem
                                         ///< (attached to some element in SubSystem)
    //Vectorize beads so the coordinates are all available in a single array.
    //@{
    static int maxbindex;//Maximum bead index alloted.
    static int vectormaxsize;//maximum number of beads that can be appended without
    // revectorization
    static int Nbeads;//Total number of beads in the system
    static vector<int> removedbindex;//stores the bead indices that have been freed


    static void revectorize(floatingpoint* coord){
        //set contiguous bindices and set coordinates.
        int idx = 0;
        for(auto b:_beads.getElements()){
            int index = 3 * idx;
            coord[index] = b->coordinate[0];
            coord[index + 1] = b->coordinate[1];
            coord[index + 2] = b->coordinate[2];
            b->_dbIndex = idx;
            idx++;
        }
        Nbeads =_beads.getElements().size();
        maxbindex = _beads.getElements().size();
        removedbindex.clear();
    }

    static void appendrevectorizeifneeded(){

        int newsize = vectormaxsize;
        //if the maximum bead index is very close to the vector size
        if(vectormaxsize - maxbindex <= bead_cache/10 )
            //new size will be increased by bead_cache
            newsize = vectormaxsize + bead_cache;
        //set parameters and revectorize
        if(newsize != vectormaxsize){
            floatingpoint *coord = CUDAcommon::serlvars.coord;
            delete[] coord;
            floatingpoint *newcoord = new floatingpoint[3 * newsize];
            CUDAcommon::serlvars.coord = newcoord;
            appendrevectorize(newcoord);
            //copyvector(newcoord, coord);
            vectormaxsize = newsize;
            //cylinder structure needs to be revecotrized as well.
            triggercylindervectorization = true;
        }
    }

	static void appendrevectorize(floatingpoint* coord){
		//set coords based on bindices.
		maxbindex = 0;
		for(auto b:_beads.getElements()){
		    maxbindex = max<int>(maxbindex, b->_dbIndex);
			int index = 3 * b->_dbIndex;
			coord[index] = b->coordinate[0];
			coord[index + 1] = b->coordinate[1];
			coord[index + 2] = b->coordinate[2];
		}
		maxbindex++;
		Nbeads =_beads.getElements().size();
	}

    //deprecated
    static void copyvector(floatingpoint* newcoord, floatingpoint* coord){
        int idx = 0;
        for(auto b:_beads.getElements()){
            idx++;
        }
    }
    //copy coodinates of this bead to the appropriate spot in coord vector.
    void  copycoordinatestovector() {
            CUDAcommon::serlvars.coord[3 * _dbIndex] = coordinate[0];
            CUDAcommon::serlvars.coord[3 * _dbIndex + 1] = coordinate[1];
            CUDAcommon::serlvars.coord[3 * _dbIndex + 2] = coordinate[2];
    }
    void resetcoordinates() {
        CUDAcommon::serlvars.coord[3 * _dbIndex] = -1.0;
        CUDAcommon::serlvars.coord[3 * _dbIndex + 1] = -1.0;
        CUDAcommon::serlvars.coord[3 * _dbIndex + 2] = -1.0;
    }
    //@}
};


#endif
