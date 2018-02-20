#ifndef MEDYAN_MembraneHierarchy_h
#define MEDYAN_MembraneHierarchy_h

#include <string>

#include "Component.h"
#include "Composite.h"

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
The closed membranes are not allowed to intersect themselves or each other.
Under this assumption, the orientable closed membranes can have a clear
hierarchy of containing relationship. The 3d space can be divided into 2
regions per membrane regardless of genus, and this hierarchy can clearly
indicate the boundaries of such regions.

This header provides a class for containing hierarchy of closed membranes
represented by a tree structure. The parent node has the membrane directly
containing the membranes of the children; the membranes of the nodes of the
same level does not contain each other.
******************************************************************************/

class MembraneHierarchy: public Composite {

private:

    Membrane* _membrane = nullptr; ///< pointer to the membrane of this level

    void printTree(string indent, bool last)const; ///< helper function for printSelf

    /**************************************************************************
    Static members
    **************************************************************************/
    static MembraneHierarchy _root; ///< The root hierarchy,
                                    ///< which should not point to any membrane

public:

    /**************************************************************************
    Ctors and Dtors
    **************************************************************************/
    MembraneHierarchy(Membrane* m): _membrane(m) {}

    /**************************************************************************
    Getters and Setters
    **************************************************************************/
    Membrane* getMembrane()const { return _membrane; }

    /**************************************************************************
    Getters and Setters (static)
    **************************************************************************/
    static MembraneHierarchy* getRoot() { return _root; }

    /**************************************************************************
    Implements Component
    **************************************************************************/
    virtual int getType()override { return 0; }
    virtual void printSelf()const override;

    /**************************************************************************
    Operations on a tree structure
    **************************************************************************/
    // When new membrane is inserted.
    // This function requires that geometry of the membrane has been updated.
    static void addMembrane(Membrane* m, MembraneHierarchy& root);
    static void addMembrane(Membrane* m) { addMembrane(m, _root); }

    // When a membrane is removed. Must be a closed membrane.
    // Returns whether something is deleted.
    static bool removeMembrane(Membrane* m, MembraneHierarchy& root);
    static bool removeMembrane(Membrane* m) { return removeMembrane(m, _root); }

};





#endif
