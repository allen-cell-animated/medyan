
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

#ifndef MEDYAN_Bubble_h
#define MEDYAN_Bubble_h

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "Composite.h"

//FORWARD DECLARATIONS
class SubSystem;
class Bead;

/// Represents a dummy point potential that is involved in mechanical equilibration. This
/// object has no chemical reactions or properties associated with it.

/*!
 *   Bubbles are artificial objects in minimization to physically represent hard spheres.
 *   They contain coordinates as well as force constants for any number of potentials. A 
 *   bubble contains a bead which is the center of the bubble, as well as a radius representing
 *   the physical size of the bubble.
 */

class Bubble : public Composite, public Trackable, public Movable, public DynamicNeighbor {

private:
    Bead* _bead;    ///< The bead representing the center of the bubble
    SubSystem* _ps; ///< The subsystem this bubble is in
    
    short _type;     ///< The type of bubble
    
    floatingpoint _radius;       ///< The radius of this bubble
    floatingpoint _kRepuls;      ///< Repulsion constant for bubble-bubble and bubble-cylinder interactions
    floatingpoint _screenLength; ///< Screening length for a repulsive potential
	floatingpoint _MTOCBendingK; ///< use for MTOC-MT bending force field
    
    int _ID;        ///< Identifier
    
    static Database<Bubble*> _bubbles; ///< Collection in SubSystem
    
    bool _isMTOC = false;   ///< If representing a MTOC
    
public:
    vector<floatingpoint> coordinate; ///< Current coordinates of bubble,
                               ///< Updated with updatePosition()
    
    /// Main constructor, sets up bead and other properties
    Bubble(SubSystem* ps, vector<floatingpoint> coordinates, short type);

    //@{
    /// Getters
    floatingpoint getRadius() {return _radius;}
    floatingpoint getRepulsionConst() {return _kRepuls;}
    floatingpoint getScreeningLength() {return _screenLength;}
	floatingpoint getMTOCBendingK() {return _MTOCBendingK;}
    
    Bead* getBead() {return _bead;}
    
    int getID() {return _ID;}
    virtual int getType() {return _type;}
    //@}
    
    void setAsMTOC() {_isMTOC = true;}
    bool isMTOC() {return _isMTOC;}

    /// Print bubble information
    virtual void printSelf();
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _bubbles.addElement(this);}
    virtual void removeFromSubSystem() {_bubbles.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Bubble*>& getBubbles() {
        return _bubbles.getElements();
    }
    /// Get the number of cylinders in this system
    static int numBubbles() {
        return _bubbles.countElements();
    }
    
    ///Update bubble position
    virtual void updatePosition();

};

#endif
