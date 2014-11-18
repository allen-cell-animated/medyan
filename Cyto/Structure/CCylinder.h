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

///FORWARD DECLARATIONS
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
    vector<unique_ptr<CMonomer>> _monomers; ///< list of monomers in this ccylinder
    
    ///REACTION CONTAINERS
    unordered_set<ReactionBase*> _internalReactions;///< list of internal reactions associated with ccylinder
    unordered_set<CCylinder*> _reactingCylinders; ///< vector of ccylinders that this ccylinder has reactions with, but not ownership
    unordered_map<CCylinder*, unordered_set<ReactionBase*>> _crossCylinderReactions; ///< map of cross-cylinder reactions owned by this ccylinder
    
    Compartment* _compartment; ///< compartment this ccylinder is in
    Cylinder* _pCylinder; ///< parent cylinder
    
    short _size = SystemParameters::Geometry().cylinderSize / SystemParameters::Geometry().monomerSize; ///<max length of full cylinder
    
public:
    ///Default constructor, sets compartment, init reactions and monomers
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
    void addCMonomer(CMonomer* monomer) { _monomers.emplace_back(unique_ptr<CMonomer>(monomer));}
    ///Get monomer at an index
    ///@note no check on index
    CMonomer* getCMonomer(int index) {return _monomers[index].get();}
    
    ///Get list of reactions associated with this CCylinder
    const unordered_set<ReactionBase*>& getInternalReactions() {return _internalReactions;}
    const unordered_set<CCylinder*>& getReactingCylinders() {return _reactingCylinders;}
    ///Get map of reactions associated with this CCylinder and others
    unordered_map<CCylinder*, unordered_set<ReactionBase*>>& getCrossCylinderReactions() {return _crossCylinderReactions;}
    
    ///REACTION MANAGEMENT FUNCTIONS
    void addInternalReaction(ReactionBase* r);
    void addCrossCylinderReaction(CCylinder* other, ReactionBase* r);
    void addReactingCylinder(CCylinder* other);

    void removeInternalReaction(ReactionBase* r);
    void removeAllInternalReactions();
    
    void removeCrossCylinderReactions(CCylinder* other);
    void removeAllCrossCylinderReactions();
    
    void removeReactingCylinder(CCylinder* other);
    void removeAllReactingCylinders();
    
//    ///Activate all reactions associated with this CCylinder
//    void activateReactions();
    
    ///Print CCylinder
    virtual void printCCylinder();
    
    ///get size of this ccylinder in number of monomers
    short getSize() {return _size;}
    
};


#endif /* defined(__CytoSim__CCylinder__) */
