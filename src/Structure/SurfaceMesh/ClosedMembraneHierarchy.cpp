#include "ClosedMembraneHierarchy.h"

#include <algorithm>
#include <stdexcept>

#include "common.h"

#include "Membrane.h"

// FORWARD DECLARATIONS
class Membrane;

void ClosedMembraneHierarchy::printSelf()const {
    cout << endl;

    cout << "ClosedMembraneHierarchy: ptr = " << this << endl;

    cout << endl;
}

void ClosedMembraneHierarchy::addMembrane(Membrane* m, const ClosedMembraneHierarchy& root) {
    auto& rootChildren = root.children();

    // Pick a point on the membrane and check the containing relationship with other membranes
    array<double, 3> p = vector2Array<double, 3>(m->getVertices().at(0)->coordinate);

    // Recursively search
    for(auto& childPtr: rootChildren) {
        
        ClosedMembraneHierarchy* hiePtr = static_cast<ClosedMembraneHierarchy*>(childPtr.get());

        if(hiePtr->_membrane == nullptr) {
            hiePtr->printSelf();
            throw runtime_error("The child node does not point to a specific membrane.");
        }

        if(hiePtr->_membrane->signedDistance(p, true) < 0) { // Is inside one of the child nodes
            ClosedMembraneHierarchy::addMembrane(m, *hiePtr); // Search that child
            return;
        }
    }

    // Now, the new membrane is outside of every child membrane.

    // First create the node under the current root (parent)
    ClosedMembraneHierarchy* newNode = new ClosedMembraneHierarchy(m);

    // Then check whether any children is inside this membrane
    for(auto& childPtr: rootChildren) {

        ClosedMembraneHierarchy* hiePtr = static_cast<ClosedMembraneHierarchy*>(childPtr.get());

        array<double, 3> hieP = vector2Array<double, 3>(hiePtr->_membrane->getVertices().at(0)->coordinate);

        if(m->signedDistance(hieP, true) < 0) { // The child membrane is inside new membrane

            // Add child to the new node
            newNode->addChild(childPtr); // Now the content of childPtr is moved and becomes the new child of the new node.
                                         // The parent of the child has also been changed.
                                         // The childPtr now should be nullptr.
            
            // For safety, check the current pointer. Don't think it is needed by standard.
            if(childPtr) {
                delete newNode; // To prevent memory leak
                throw logic_error("Original child pointer not null after being moved.");
            }
            
        }
    }

    // Then remove all null children from root using erase-remove idiom.
    // Good news: unique_ptr can be compared to nullptr_t.
    rootChildren.erase(remove(rootChildren.begin(), rootChildren.end(), nullptr), rootChildren.end());

    // Finally add the new node to the current root
    root.addChild(newNode); // Also manages the deletion of newNode

}
