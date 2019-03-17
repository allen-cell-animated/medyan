
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

class Bubble : public Composite, public Trackable, public Movable, public DynamicNeighbor,
    public Database< Bubble > {

private:
    Bead* _bead;    ///< The bead representing the center of the bubble
    SubSystem* _ps; ///< The subsystem this bubble is in
    
    short _type;     ///< The type of bubble
    
    double _radius;       ///< The radius of this bubble
    double _kRepuls;      ///< Repulsion constant for bubble-bubble and bubble-cylinder interactions
    double _screenLength; ///< Screening length for a repulsive potential
    double _MTOCBendingK; ///< use for MTOC-MT bending force field
    
    bool _isMTOC = false;   ///< If representing a MTOC
    
public:
    vector<double> coordinate; ///< Current coordinates of bubble,
                               ///< Updated with updatePosition()
    
    /// Main constructor, sets up bead and other properties
    Bubble(SubSystem* ps, vector<double> coordinates, short type);

    //@{
    /// Getters
    double getRadius() {return _radius;}
    double getRepulsionConst() {return _kRepuls;}
    double getScreeningLength() {return _screenLength;}
    double getMTOCBendingK() {return _MTOCBendingK;}
    
    Bead* getBead() {return _bead;}
    
    virtual int getType() {return _type;}
    //@}
    
    void setAsMTOC() {_isMTOC = true;}
    bool isMTOC() {return _isMTOC;}

    /// Print bubble information
    virtual void printSelf();
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { }
    virtual void removeFromSubSystem() {}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Bubble*>& getBubbles() {
        return getElements();
    }
    /// Get the number of cylinders in this system
    static int numBubbles() {
        return getElements().size();
    }
    
    ///Update bubble position
    virtual void updatePosition();

};

#endif
