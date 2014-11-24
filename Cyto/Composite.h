//
//  Composite.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/29/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Composite_h
#define CytoSim_Composite_h

#include <string>
#include <vector>
#include <deque>
#include <typeinfo>

#include "common.h"

#include "Component.h"
#include "Species.h"

///FORWARD DECLARATIONS
class Visitor;
class SpeciesVisitor;
class ReactionVisitor;


/// Composite class is the aggregating class for the Composite pattern

/*! The Composite pattern allows building of complex hieararchical objects, with convenient
 *  methods for applying a function to all nodes (i.e. the Visitor pattern). Each node in the
 *  hieararchy may have a parent and may contain several children nodes. A class that is derived
 *  from Composite can contain children Component objects.
 *  @note Composite objects may contain Species and ReactionBase collections, however, this is
 *  treated seperately from the Composite pattern (i.e. separate methods exist for the corresponding
 *  access to elements, changes, etc.)
 */    
class Composite : public Component {
private:
    vector<unique_ptr<Component>> _children;///< Child node pointers of this Composite node

public:
    /// Default Constructor does nothing.
    Composite() :  Component() {}
    
    /// Virtual Destructor does nothing.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~Composite() noexcept {}
    
    /// Implements the apply_if() method of the Component class by recursively applying it
    /// to itself and all its children.
    virtual bool apply (Visitor &v) override;
    
    /// Implements the apply_if() method of the Component class by recursively applying it
    /// to itself and all its children that contain Species.
    virtual bool apply (SpeciesVisitor &v) override;
    
    /// Implements the apply_if() method of the Component class by recursively applying it
    /// to itself and all its children that contain ReactionBase.
    virtual bool apply (ReactionVisitor &v) override;
    
    /// Returns true.
    virtual bool isComposite() const override {return true;}
    
    /// Returns the full name of this node.
    virtual string getFullName() const override {return "Composite";}; 
    
    /// Adds a Component child to this Composite node
    /// @param child - is a unique_ptr<Component>, hence, this node takes the memory ownership
    /// of the corresponding child pointer.
    /// @note the rvalue semantics - the unique_ptr<...> cannot be copied, but only can be moved
    virtual void addChild (unique_ptr<Component> &&child) {
        _children.push_back(move(child));
        _children.back()->setParent(this);
    }
    
    
    /// Remove *child from this node's children. The child's destructor will be called and the memory will be freed.
    virtual void removeChild (Component* child) {
        auto child_iter = find_if(_children.begin(),_children.end(),
                    [&child](const unique_ptr<Component> &element)
                     {
                         return element.get()==child ? true : false;
                     });
        if(child_iter!=_children.end())
            _children.erase(child_iter);
        else
            throw out_of_range("Composite::removeChild(): The name child not found");
    }
    
    
    /// Returns the number of immediate children of this node.
    /// @note Species and reactions and not included in this tally
    virtual size_t numberOfChildren() const {return children().size();}

    /// Returns the number of Species being immediately managed by this node (i.e. not counting Species
    /// belonging to children nodes, etc.
    virtual size_t numberOfSpecies () const override {return 0;}

    /// Return the total number of nodes contained under this node's hieararchy
    /// @note This is a recursive call, and all nodes under this node are visited.
    virtual size_t countDescendents() const override {
        size_t sum = numberOfChildren();
        for (auto &c : children())
            sum+=c->numberOfChildren();
        return sum;
    }

    
    /// Returns the number of Species being managed by this node and its descendent nodes
    virtual size_t countSpecies() const override {
        size_t sum = 0;
        if(this->isSpeciesContainer())
            sum+=this->numberOfSpecies();
        for (auto &c : children())
            sum+=c->countSpecies();
        return sum;
    }
    
    /// Returns the number of ReactionBase objets being immediately managed by this node (i.e. not counting reactions
    /// belonging to children nodes, etc.
    virtual size_t numberOfReactions() const override {return 0;}
    
    /// Returns the number of ReactionBase objects being managed by this node and its descendent nodes
    virtual size_t countReactions() const override {
        size_t sum = 0;
        if(this->isReactionsContainer())
            sum+=this->numberOfReactions();
        for (auto &c : children())
            sum+=c->countReactions();
        return sum;
    }
    
    /// Returns a reference to the container of Component children of this node
    virtual vector<unique_ptr<Component>>& children () {return _children;}

    /// Returns a const reference to the container of Component children of this node
    virtual const vector<unique_ptr<Component>>& children () const {return _children;}

    /// Returns a pointer to the i-th Component child of this node
    virtual Component* children (size_t i) {return _children[i].get();}
    
};
       
#endif
