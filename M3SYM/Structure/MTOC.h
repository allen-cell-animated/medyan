
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

#ifndef M3SYM_MTOC_h
#define M3SYM_MTOC_h

#include "common.h"

#include "Database.h"
#include "Trackable.h"

//FORWARD DECLARATIONS
class Bubble;
class Filament;

///A class to represent the structure of a microtubule organizing center (MTOC)
/*!
 *  The MTOC class provides a base structure for a mechanical MTOC, which includes
 *  a Bubble representing the MTOC and [Filaments](@ref Filament) which are 
 *  attached to it at the minus end by some sort of potential.
 *
 *  This class only has functions to set and get the constituent Bubble and Filament
 *  and is only mechanically relevant for now, but may be extended to have chemical 
 *  properties in the future.
 */
class MTOC : public Trackable {
    
private:
    Bubble* _bubble; ///< A bubble that physically represents the MTOC
    vector<Filament*> _filaments; ///< An ordered vector of filaments in the MTOC
    
    static Database<MTOC*> _mtocs; ///< Collection of beads in SubSystem
public:
    ///Constructor
    MTOC() : Trackable() {}
    
    //@{
    ///Setters
    void setBubble(Bubble* b) {_bubble = b;}
    
    void addFilament(Filament* f) {_filaments.push_back(f);}
    //@}
    
    //@{
    ///Getters
    Bubble* getBubble() {return _bubble;}
    const vector<Filament*>& getFilaments() {return _filaments;}
    //@}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _mtocs.addElement(this);}
    virtual void removeFromSubSystem() {_mtocs.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<MTOC*>& getMTOCs() {
        return _mtocs.getElements();
    }
    /// Get the number of MTOCs in this system
    static int numMTOCs() {
        return _mtocs.countElements();
    }
    
};

#endif
