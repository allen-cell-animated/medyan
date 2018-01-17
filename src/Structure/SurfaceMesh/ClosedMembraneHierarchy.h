#ifndef MEDYAN_ClosedMembraneHierarchy_h
#define MEDYAN_ClosedMembraneHierarchy_h

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

class ClosedMembraneHierarchy: public Composite {

private:

    Membrane* _membrane = nullptr; // pointer to the membrane of this level

public:

    /**************************************************************************
    Ctors and Dtors
    **************************************************************************/
    ClosedMembraneHierarchy(Membrane* m): _membrane(m) {}

    /**************************************************************************
    Implements Component
    **************************************************************************/
    virtual int getType()override { return 0; }
    virtual void printSelf()const override;

    /**************************************************************************
    Operations on a tree structure
    **************************************************************************/
    // When new membrane is inserted. Must be a closed membrane.
    static void addMembrane(Membrane* m, const ClosedMembraneHierarchy& root);

    // When a membrane is removed. Must be a closed membrane.
    static void removeMembrane(Membrane* m, const ClosedMembraneHierarchy& root);

};





#endif
