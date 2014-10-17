//
//  CCylinder.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CCylinder__
#define __CytoSim__CCylinder__

#include <iostream>

#include "common.h"
#include "Compartment.h"
#include "CMonomer.h"
#include "SystemParameters.h"

class Cylinder;

/// CCylinder class holds all CMonomers and reactions associated with its species
/*! 
 *  The CCylinder class has lists of the monomers and bounds that it contains.
 *
 *  A CCylinder also stores the reaction associated with it internally, as well as cross-cylinder
 *  reactions (still "owned" by this ccylinder, but has reactants and/or products that contain a
 *  species in another ccylinder). This information is stored in the reaction map. This map can be
 *  used to add and delete reactions when a ccylinder changes compartment or is removed from the
 *  system.
 *
 *  Lastly, the CCylinder stores a list of other ccylinders that it has cross-cylinder reactions with,
 *  but does not have ownership of these reactions. This list can be used to change or delete reactions
 *  when needed, but will not be tied to a compartment transfer of this ccylinder.
 */
class CCylinder {
    
protected:
    std::vector<std::unique_ptr<CMonomer>> _monomers; ///< list of monomers in this ccylinder
    
    ///REACTION CONTAINERS
    std::unordered_set<ReactionBase*> _internalReactions;///< list of internal reactions associated with ccylinder
    std::unordered_set<CCylinder*> _reactingCylinders; ///< vector of ccylinders that this ccylinder has reactions with, but not ownership
    std::unordered_map<CCylinder*, std::vector<ReactionBase*>> _crossCylinderReactions; ///< map of cross-cylinder reactions owned by this ccylinder
    
    Compartment* _compartment; ///< compartment this ccylinder is in

    Cylinder* _pCylinder; ///< parent cylinder
    
    short _size = SystemParameters::Geometry().cylinderSize / SystemParameters::Geometry().monomerSize; ///<max length of full cylinder

public:
    ///Default constructor, sets compartment, init reaction containers
    CCylinder(Compartment* c) : _compartment(c) {}
    
    /// Copy constructor
    /// @note This constructor will create a new ccylinder with different species and reactions
    /// within the compartment that is chosen as a parameter to the constructor. The copied and
    /// original ccylinders will not share reactions or species, but will be copied into a new compartment.
    CCylinder(const CCylinder& rhs, Compartment* c);
    
    /// Assignment is not allowed
    CCylinder& operator=(CCylinder &rhs) = delete;
    
    ///Default destructor, explicitly removes monomers and bounds (including their species, rxns)
    /// Removes all reactions associated with this ccylinder, including ones owned by this as well
    /// as other ccylinders.
    ~CCylinder();
    
    ///Clone, calls copy constructor
    virtual CCylinder* clone(Compartment* c) {
        return new CCylinder(*this, c);
    }
    
    ///get filament compartment
    Compartment* getCompartment() {return _compartment;}
    
    ///set parent cylinder
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    ///get parent cylinder
    Cylinder* getCylinder() {return _pCylinder;}
    
    ///Add a monomer to this CCylinder
    void addCMonomer(CMonomer* monomer) { _monomers.emplace_back(std::unique_ptr<CMonomer>(monomer));}
    ///Get monomer at an index
    ///@note no check on index
    CMonomer* getCMonomer(int index) {return _monomers[index].get();}
    
    ///Add an internal reaction to this CCylinder
    void addInternalReaction(ReactionBase* r);
    ///Add a reaction across this ccylinder and another ccylinder
    ///Adds to this reaction map as well as the other ccylinder's reacting cylinder list
    void addCrossCylinderReaction(CCylinder* other, ReactionBase* r);
    ///Add a reacting ccylinder to this one
    void addReactingCylinder(CCylinder* other);
    
    ///remove an internal reaction from this ccylinder
    ///@note no check on whether r is in the reactions list
    void removeInternalReaction(ReactionBase* r);
    ///remove all internal reactions from this ccylinder
    void removeAllInternalReactions();
    
    ///clear all reactions involving another ccylinder
    ///@note also removes reaction from other ccylinder's reacting cylinder lsist
    void removeCrossCylinderReactions(CCylinder* other);
    ///clear all reactions involving other cylinders
    ///@note also removes all reactions from other ccylinder's reacting cylinder list
    void removeAllCrossCylinderReactions();
    
    ///remove a reacting ccylinder from this one,
    ///@note - also removes the other's cross-cylinder reactions involving this ccylinder
    void removeReactingCylinder(CCylinder* other);
    ///remove all reacting ccylinders reactions,
    ///@note - similar to above but for all in reacting cylinder list
    void removeAllReactingCylinders();
    
    ///Get list of reactions associated with this CCylinder
    const std::unordered_set<ReactionBase*>& getInternalReactions() {return _internalReactions;}
    const std::unordered_set<CCylinder*>& getReactingCylinders() {return _reactingCylinders;}
    ///Get map of reactions associated with this CCylinder and others
    std::unordered_map<CCylinder*, std::vector<ReactionBase*>>& getCrossCylinderReactions() {return _crossCylinderReactions;}
    
    ///Update all reactions associated with this CCylinder
    void updateReactions();
    ///Print CCylinder
    virtual void printCCylinder();
    
    ///get size of this ccylinder in number of monomers
    short size() {return _size;}
    
};


#endif /* defined(__CytoSim__CCylinder__) */
