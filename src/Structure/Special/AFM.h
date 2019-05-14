
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

#ifndef MEDYAN_AFM_h
#define MEDYAN_AFM_h

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Composite.h"

//FORWARD DECLARATIONS
class Bubble;
class Filament;

///A class to represent the structure of a AFM attachemnt system
/*!
 *  The AFM class provides a base structure for a mechanical AFM, which includes
 *  a Bubble (or two bubbles) representing the AFM and [Filaments](@ref Filament) which are
 *  attached to it at the minus end by some sort of potential.
 *
 *  This class only has functions to set and get the constituent Bubble and Filament
 *  and is only mechanically relevant for now, but may be extended to have chemical
 *  properties in the future.
 */
class AFM : public Composite, public Trackable {
    
private:
    Bubble* _bubble; ///< A bubble that physically represents the AFM bubble
    vector<Filament*> _filaments; ///< An ordered vector of filaments in the AFM bubble
    
    int _ID; ///< Unique identifier
    
    static Database<AFM*> _afms; ///< Collection of beads in SubSystem
public:
    ///Constructor
    AFM() : Trackable(), _ID(_afms.getID()) {}
    
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
    virtual void addToSubSystem() { _afms.addElement(this);}
    virtual void removeFromSubSystem() {_afms.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<AFM*>& getAFMs() {
        return _afms.getElements();
    }
    /// Get the number of AFMs in this system
    static int numAFMs() {
        return _afms.countElements();
    }
    
    virtual void printSelf();
    
    //GetType implementation just returns zero (no AFM types yet)
    virtual int getType() {return 0;}
    
};

#endif

