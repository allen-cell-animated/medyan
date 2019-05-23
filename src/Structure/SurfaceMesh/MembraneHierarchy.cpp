#include "Structure/SurfaceMesh/MembraneHierarchy.h"

#include <algorithm>
#include <stdexcept>

#include "common.h"
#include "MathFunctions.h"
using namespace mathfunc;

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Vertex.h"

// Forward declarations
class Membrane;

// Static member initialization
MembraneHierarchy MembraneHierarchy::_root(nullptr);

// Member function definitions
void MembraneHierarchy::printSelf()const {
    cout << endl;

    cout << "MembraneHierarchy: ptr = " << this << endl;

    cout << endl << "Tree structure:" << endl;
    printTree("", true);

    cout << endl;
}

void MembraneHierarchy::printTree(string indent, bool last)const {
    cout << indent;
    if (last) {
        cout << "\\-";
        indent += "  ";
    }
    else {
        cout << "|-";
        indent += "| ";
    }
    
    cout << this;
    if(_membrane)
        cout << " (Mem Ptr: " << _membrane << " Id: " << _membrane->getId() << ")";
    else
        cout << " (No membrane attached)";
    
    cout << endl;

    size_t n = numberOfChildren();
    for (size_t idx = 0; idx < n; ++idx)
        static_cast<const MembraneHierarchy*>(children(idx))->printTree(indent, idx == n - 1);
}

void MembraneHierarchy::addMembrane(Membrane* m, MembraneHierarchy& root) {
    auto& rootChildren = root.children();

    // Pick a point on the membrane and check the containing relationship with other membranes
    const mathfunc::Vec3 p = m->getMesh().getVertexAttribute(0).getCoordinate();

    // Recursively search
    for(auto& childPtr: rootChildren) {
        
        MembraneHierarchy* hiePtr = static_cast<MembraneHierarchy*>(childPtr.get());

        if(hiePtr->_membrane == nullptr) {
            hiePtr->printSelf();
            throw runtime_error("The child node does not point to a specific membrane.");
        }

        if(hiePtr->_membrane->isClosed() && hiePtr->_membrane->contains(p)) { // Is inside one of the child nodes
            MembraneHierarchy::addMembrane(m, *hiePtr); // Search that child
            return;
        }
    }

    // Now, the new membrane is outside of every child membrane.

    // First create a new node
    MembraneHierarchy* newNode = new MembraneHierarchy(m);

    // Then check whether any children is inside this membrane.
    // Open membranes cannot contain any children.
    if(m->isClosed()) {

        for(auto& childPtr: rootChildren) {

            MembraneHierarchy* hiePtr = static_cast<MembraneHierarchy*>(childPtr.get());

            const mathfunc::Vec3 hieP = hiePtr->_membrane->getMesh().getVertexAttribute(0).getCoordinate();

            if(m->contains(hieP)) { // The child membrane is inside new membrane

                // Add child to the new node
                newNode->addChild(move(childPtr)); // Now the content of childPtr is moved and becomes the new child of the new node.
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
    }

    // Finally add the new node to the current root
    root.addChild(unique_ptr<Component>(newNode)); // Also manages the deletion of newNode

}

bool MembraneHierarchy::removeMembrane(Membrane* m, MembraneHierarchy& root) {
    // Find the membrane to be removed by recursive search
    // Only 1 node will be deleted

    MembraneHierarchy* nodeToBeDeleted = nullptr;
    
    for(auto& childPtr: root.children()) {
        
        MembraneHierarchy* hiePtr = static_cast<MembraneHierarchy*>(childPtr.get());

        if(hiePtr->_membrane == m) { // Found the node to be removed
            nodeToBeDeleted = hiePtr;
        }
        else {
            if(MembraneHierarchy::removeMembrane(m, *hiePtr)) return true;
            // else the search continues
        }
    }

    if(nodeToBeDeleted) {

        // Bring all its children under its parent.
        for(auto& childHiePtr: nodeToBeDeleted->children()) {
            root.addChild(move(childHiePtr));
        }

        // Remove the node from root and destroy it.
        root.removeChild(nodeToBeDeleted);

        // Exit function
        return true;

    } else {
        return false;
    }

}
