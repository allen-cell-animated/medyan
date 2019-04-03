
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

#ifndef MEDYAN_MTOC_h
#define MEDYAN_MTOC_h

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Composite.h"

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
class MTOC : public Composite, public Trackable,
    public Database< MTOC, false > {
    
private:
    Bubble* _bubble; ///< A bubble that physically represents the MTOC
    vector<Filament*> _filaments; ///< An ordered vector of filaments in the MTOC
    
public:
    ///Constructor
    MTOC() : Trackable() {}
    
    //@{
    ///Setters
    void setBubble(Bubble* b);
    
    void addFilament(Filament* f) {_filaments.push_back(f);}
    //@}
    
    //@{
    ///Getters
    Bubble* getBubble() {return _bubble;}
    const vector<Filament*>& getFilaments() {return _filaments;}
    //@}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() { }
    virtual void removeFromSubSystem() {}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<MTOC*>& getMTOCs() {
        return getElements();
    }
    /// Get the number of MTOCs in this system
    static int numMTOCs() {
        return getElements().size();
    }
    
    virtual void printSelf();
    
    //GetType implementation just returns zero (no MTOC types yet)
    virtual int getType() {return 0;}
    
};

#endif
