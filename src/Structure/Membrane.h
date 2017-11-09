#ifndef MEDYAN_Membrane_h
#define MEDYAN_Membrane_h

#include <vector>

#include "Trackable.h"
#include "Composite.h"

// FORWARD DECLARATIONS
class SubSystem;
class Triangle;
class Bead;

class Membrane: public Composite, public Trackable {

    friend class Controller;

private:

    SubSystem* _subSystem; // SubSystem pointer

    static Database<Membrane*> _membranes; // Collection in SubSystem

public:
    // SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _membranes.addElement(this); }
    virtual void removeFromSubSystem() { _membranes.removeElement(this); }
    
    /// Get all instances of this class from the SubSystem
    static const vector<Membrane*>& getMembranes() {
        return _membranes.getElements();
    }
    /// Get the number of membranes in this system
    static int numFilaments() {
        return _membranes.countElements();
    }


};








#endif
