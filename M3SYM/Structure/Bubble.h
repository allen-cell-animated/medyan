
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

#ifndef M3SYM_Bubble_h
#define M3SYM_Bubble_h

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
    
    double _radius;       ///< The radius of this bubble
    double _kRepuls;      ///< Repulsion constant for bubble-bubble and bubble-cylinder interactions
    double _screenLength; ///< Screening length for a repulsive potential
    
    int _ID;  ///< Identifier
    
    static Database<Bubble*> _bubbles; ///< Collection in SubSystem
    
public:
    vector<double> coordinate; ///< Current coordinates of bubble,
                               ///< Updated with updatePosition()
    
    /// Main constructor, sets up bead and other properties
    Bubble(SubSystem* ps, vector<double> coordinates,
           short type, double radius, double kRepuls, double screenLength);

    //@{
    /// Getters
    double getRadius() {return _radius;}
    double getRepulsionConst() {return _kRepuls;}
    double getScreeningLength() {return _screenLength;}
    
    Bead* getBead() {return _bead;}
    
    int getID() {return _ID;}
    //@}
    
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
