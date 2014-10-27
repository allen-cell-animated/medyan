//
//  BoundaryElement.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElement__
#define __Cyto__BoundaryElement__

#include <vector>
#include <iostream>

#include "common.h"
#include "Gcontroller.h"

class Bead;

///BoundaryElement class represents an element of a boundary surface
/*!
 * The BoundaryElement class is a representation of a boundary surface element, which can interact
 * with other elements in the system, including other BoundaryElements as well as Beads in filaments.
 * Together, a collection of boundary elements make up a BoundarySurface.
 */
class BoundaryElement {
    
protected:
    
    Compartment* _compartment; ///< Compartment that this boundary element is currently in
    std::vector<Bead*> _beads; ///< Beads that this boundary element could interact with
    
    std::vector<double> _coords; ///< coordinates of this boundary element
    
public:
    ///Default constructor
    BoundaryElement(std::vector<double> coords) : _coords(coords) {
    
        ///set the compartment given the initial coordinates
        try {_compartment = GController::getCompartment(coords);}
        catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    }
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~BoundaryElement() noexcept {
        ///remove from compartment
        _compartment->removeBoundaryElement(this);
    }
    
    ///add a bead to list of interacting beads
    void addBead(Bead* b) {_beads.push_back(b);}
    ///Remove a bead from list of interacting beads
    ///@note does nothing if bead is not in interacting list already
    void removeBead(Bead* b) {
        auto it = std::find(_beads.begin(), _beads.end(), b);
        if(it != _beads.end()) _beads.erase(it);
    }
    
    ///Get the compartment that this element is in
    Compartment* getCompartment() {return _compartment;}
    
    ///return coordinates of boundary element
    const std::vector<double>& coords() {return _coords;}
    ///return set of beads
    const std::vector<Bead*>& beads() {return _beads;}
    
    ///Implement for all boundary elements
    virtual double distance(const std::vector<double>& point) = 0;
    virtual double stretchedDistance(const std::vector<double>& point, const std::vector<double>& force, double d) = 0;
    virtual const std::vector<double> normal(const std::vector<double> &point) = 0;
    
    virtual double getRepulsionConst() = 0;
    virtual double getScreeningLength() = 0;
 
};



#endif /* defined(__Cyto__BoundaryElement__) */
