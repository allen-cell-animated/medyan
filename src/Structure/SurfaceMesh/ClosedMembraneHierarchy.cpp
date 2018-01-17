#include "ClosedMembraneHierarchy.h"

#include <algorithm>
#include <stdexcept>

#include "common.h"

#include "Membrane.h"
#include "Vertex.h"

// FORWARD DECLARATIONS
class Membrane;

void ClosedMembraneHierarchy::printSelf()const {
    cout << endl;

    cout << "ClosedMembraneHierarchy: ptr = " << this << endl;

    printTree("", true);

    cout << endl;
}

void ClosedMembraneHierarchy::printTree(string indent, bool last)const {
    cout << indent;
    if (last) {
        cout << "\-";
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
        children()[idx]->printTree(indent, idx == n - 1);
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

    // First create a new node
    ClosedMembraneHierarchy* newNode = new ClosedMembraneHierarchy(m);

    // Then check whether any children is inside this membrane
    for(auto& childPtr: rootChildren) {

        ClosedMembraneHierarchy* hiePtr = static_cast<ClosedMembraneHierarchy*>(childPtr.get());

        array<double, 3> hieP = vector2Array<double, 3>(hiePtr->_membrane->getVertices().at(0)->coordinate);

        if(m->signedDistance(hieP, true) < 0) { // The child membrane is inside new membrane

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

    // Finally add the new node to the current root
    root.addChild(unique_ptr<Component>(newNode)); // Also manages the deletion of newNode

}

bool ClosedMembraneHierarchy::removeMembrane(Membrane* m, const ClosedMembraneHierarchy& root) {
    // Find the membrane to be removed by recursive search
    // Only 1 node will be deleted

    ClosedMembraneHierarchy* nodeToBeDeleted = nullptr;
    
    for(auto& childPtr: root.children()) {
        
        ClosedMembraneHierarchy* hiePtr = static_cast<ClosedMembraneHierarchy*>(childPtr.get());

        if(hiePtr->_membrane == m) { // Found the node to be removed
            nodeToBeDeleted = hiePtr;
        }
        else {
            if(ClosedMembraneHierarchy::removeMembrane(m, *hiePtr)) return true;
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
