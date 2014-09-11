//
//  Visitor.h
//  CytoSim
//
//  Created by Garegin Papoian on 6/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Visitor_h
#define CytoSim_Visitor_h

#include <iostream>
#include "Component.h"

/// An abstract interface to traverse those nodes in the Composite pattern which fulfil a certain predicate.

/*! The visitor pattern allows a functor to visit each node of the Composite pattern. Concrete
 *  classes derived from Visitor must implement the visit(Component *c) method. Visitor, v, allows a
 *  predicate to be applied at each node, *c, and cv->visit(c) is called only if that node
 *  satisfies the predicate. For example, this can be used to visit only filaments in the
 *  hieararchy.
 */
class Visitor {
public:
    /// When this conditional visitor, *cv, is propagated through the Composite hieararchy, at each
    /// Component node pointer *c, the following method is called: v->visit(c), if v->pred(c) returns true.
    /// If v->visit(c) returns a false value,  then the further propagation of the tree is terminated.
    bool visit(Component *c) {
        if(predImpl(c)){
            return visitImpl(c);
        }
        else{
            return true;
        }
    }
    
    /// Virtual destructor
    virtual ~Visitor() {}
    
protected:
    /// When this visitor, *cv, is propagated through the Composite hieararchy, at each
    /// Component node pointer *c, the following method is called: cv->visit(c).
    /// If a false value is returned from this function called, then the further
    /// propagation of the tree is terminated.
    virtual bool visitImpl(Component *c) = 0;
    
    /// Return true if the Component *c satisfies the desired predicate
    virtual bool predImpl(Component *c) {return true;};
};

        
/// The VisitorLambda allows using C++11 lambda expressions to set the action to be performed on each node, and also check via a lambda predicate whether the given node needs to be acted upon.
class VisitorLambda : public Visitor {
public:
    /// Sets the action to be perfomed on each node as std::function<bool (Component*)> g
    void setLambda(std::function<bool (Component*)> g) {_f=g;}
    
    /// Sets the predicate which checks whether the node should be evaluated as std::function<bool (Component*)> g
    void setLambdaPred(std::function<bool (Component*)> pred) {_pred=pred;}
protected:
    std::function<bool (Component*)> _f;///< The function to be perfomed on each node
    std::function<bool (Component*)> _pred = [](Component*){return true;};///< The predicate which checks whether the node should be evaluated
    
    /// Implementation of the visit(...) method
    virtual bool visitImpl(Component *c) override {return _f(c);}
    
    /// Implementation of the pred(...) method
    virtual bool predImpl(Component *c) override {return _pred(c);}
};
    
template <class ComponentType>
class TypeVisitor : public Visitor {
protected:
    virtual bool predImpl(Component *c) override {
        return isSame<ComponentType>(*c);
    }
};

template <class ComponentType>
class TypeVisitorLambda : public VisitorLambda {
protected:
    virtual bool predImpl(Component *c) override {
        return isSame<ComponentType>(*c) && VisitorLambda::predImpl(c);
    }
};

/*! The visitor pattern allows a functor to visit each node of the Composite pattern. Concrete
 *  classes derived from Visitor must implement the visit(Component *c) method.
 */
class SpeciesVisitor {
public:
    /// When this visitor, *v, is propagated through the Composite hieararchy, at each
    /// Component node pointer *c, the following method is called: v->visit(s), for each of the
    /// Species *s of that node.
    /// If a false value is returned from this function called, then the further
    /// propagation of the tree is terminated.
    bool visit(Species *s)
    {
        if(predImpl(s)){
            return visitImpl(s);
        }
        else{
            return true;
        }
        
    }

    /// Virtual destructor
    virtual ~SpeciesVisitor() {}
protected:
    /// When this visitor, *cv, is propagated through the Composite hieararchy, at each
    /// Component node pointer *c, the following method is called: cv->visit(c).
    /// If a false value is returned from this function called, then the further
    /// propagation of the tree is terminated.
    virtual bool visitImpl(Species *s) = 0;
    
    /// Return true if the Species *s satisfies the desired predicate
    virtual bool predImpl(Species *s) {return true;}
};



/// The SpeciesVisitorLambda allows using C++11 lambda expressions to set the action to be performed on
/// each Species of the node and its descendent nodes
class SpeciesVisitorLambda : public SpeciesVisitor {
public:
    /// Sets the action to be perfomed on each node as std::function<bool (Component*)> g
    void setLambda(std::function<bool (Species *)> g) {_f=g;}
    
    /// Sets the predicate which checks whether the node should be evaluated as std::function<bool (Component*)> g
    void setLambdaPred(std::function<bool (Species*)> pred) {_pred=pred;}
protected:
    std::function<bool (Species *)> _f;///< The function to be perfomed on each node
    std::function<bool (Species *)> _pred =[](Species *){return true;};///< The predicate which checks whether the Species should be evaluated
    
    /// Implementation of the visit(...) method
    virtual bool visitImpl(Species *s) override {return _f(s);}
    
    /// Implementation of the pred(...) method
    virtual bool predImpl(Species *s) override {return _pred(s);}
};


/*! The visitor pattern allows a functor to visit each node of the Composite pattern. Concrete
 *  classes derived from Visitor must implement the visit(Component *c) method.
 */
class ReactionVisitor {
public:
    /// When this visitor, *v, is propagated through the Composite hieararchy, at each
    /// Component node pointer *c, the following method is called: v->visit(s), for each of the
    /// Reaction *r of that node.
    /// If a false value is returned from this function called, then the further
    /// propagation of the tree is terminated.
    bool visit(ReactionBase *r)
    {
        //        std::cout << "ConditionalVisitor::visit_if: checking " << c->getFullName() << std::endl;
        if(predImpl(r)){
            return visitImpl(r);
        }
        else{
            //            std::cout << "ConditionalVisitor::visit_if: Predicate failed, moving on..." << std::endl;
            return true;
        }
        
    }
    
    /// Virtual destructor
    virtual ~ReactionVisitor() {}
protected:
    /// When this visitor, *cv, is propagated through the Composite hieararchy, at each
    /// Component node pointer *c, the following method is called: cv->visit(c).
    /// If a false value is returned from this function called, then the further
    /// propagation of the tree is terminated.
    virtual bool visitImpl(ReactionBase *r) = 0;
    
    /// Return true if the Reaction *r satisfies the desired predicate
    virtual bool predImpl(ReactionBase *r) = 0;
};



/// The ReactionVisitorLambda allows using C++11 lambda expressions to set the action to be performed on
/// each Reaction of the node and its descendent nodes
class ReactionVisitorLambda : public ReactionVisitor {
public:
    /// Default Constructor
    ReactionVisitorLambda() = default;
    
    /// Constructor:
    /// @param std::function<bool (Reaction*)> g defines the action to be perfomed on each Reaction of the node
    /// @param std::function<bool (Reaction*)> pred defines the predicate whether the particular Reaction should be processed
    ReactionVisitorLambda(std::function<bool (ReactionBase *)> g, std::function<bool (ReactionBase *)> pred = [](ReactionBase *){return true;}) : ReactionVisitor(), _f(g), _pred(pred) {}
    
    /// Sets the action to be perfomed on each node as std::function<bool (Component*)> g
    void setLambda(std::function<bool (ReactionBase *)> g) {_f=g;}
    
    /// Sets the predicate which checks whether the node should be evaluated as std::function<bool (Component*)> g
    void setLambdaPred(std::function<bool (ReactionBase*)> pred) {_pred=pred;}
protected:
    std::function<bool (ReactionBase *)> _f;///< The function to be perfomed on each node
    std::function<bool (ReactionBase *)> _pred;///< The predicate which checks whether the Reaction should be evaluated
    
    /// Implementation of the visit(...) method
    virtual bool visitImpl(ReactionBase *r) override {return _f(r);}
    
    /// Implementation of the pred(...) method
    virtual bool predImpl(ReactionBase *r) override {return _pred(r);}
};

#endif
