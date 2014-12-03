
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_CCylinder_h
#define M3SYM_CCylinder_h

#include "common.h"

#include "CMonomer.h"
#include "Compartment.h"

#include "SystemParameters.h"

//FORWARD DECLARATIONS
class Cylinder;

/// Holds all [CMonomers](@ref CMonomer) and [Reactions](@ref Reaction) associated with it.
/*! 
 *  The CCylinder class has lists of the [CMonomers](@ref CMonomer) that it contains.
 *
 *  A CCylinder also stores the [Reactions] (@ref Reaction) associated with it internally, as well as cross-cylinder
 *  reactions (still "owned" by this CCylinder, but has reactants and/or products that contain a
 *  species in another CCylinder). This information is stored in the reaction map. This map can be
 *  used to add and delete reactions when a CCylinder changes compartment or is removed from the
 *  system.
 *
 *  Lastly, the CCylinder stores a list of other CCylinders that it has cross-cylinder reactions with,
 *  but does not have ownership of these reactions. This list can be used to change or delete reactions
 *  when needed, but will not be tied to a Compartment transfer of this CCylinder.
 */
class CCylinder {
    
private:
    vector<unique_ptr<CMonomer>> _monomers; ///< List of monomers
    
    ///REACTION CONTAINERS
    unordered_set<ReactionBase*> _internalReactions;///< List of internal reactions associated
    unordered_set<CCylinder*> _reactingCylinders; ///< Set of ccylinders that this ccylinder has reactions with, but not ownership
    unordered_map<CCylinder*, unordered_set<ReactionBase*>> _crossCylinderReactions; ///< Map of cross-cylinder reactions owned
    
    Compartment* _compartment; ///< Compartment this ccylinder is in
    Cylinder* _pCylinder; ///< Parent cylinder
    
    short _size = SystemParameters::Geometry().cylinderSize / SystemParameters::Geometry().monomerSize; ///< Max length of full Cylinder
    
public:
    /// Default constructor, sets compartment
    CCylinder(Compartment* c) : _compartment(c) {}
    
    /// Copy constructor
    /// @note This constructor will create a new CCylinder with different Species and [Reactions](@ref Reaction)
    /// within the Compartment that is chosen as a parameter to the constructor. The copied and
    /// original CCylinder will not share reactions or species, but will be copied into a new Compartment.
    CCylinder(const CCylinder& rhs, Compartment* c);
    
    /// Assignment is not allowed
    CCylinder& operator=(CCylinder &rhs) = delete;
    
    /// Default destructor, explicitly removes CMonomers(including their species, rxns)
    /// Removes all reactions associated with this CCylinder, including ones owned by this as well
    /// as other CCylinders.
    ~CCylinder();
    
    /// Clone, calls copy constructor
    virtual CCylinder* clone(Compartment* c) {
        return new CCylinder(*this, c);
    }
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Set parent
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    /// Get parent
    Cylinder* getCylinder() {return _pCylinder;}
    
    /// Add a monomer
    void addCMonomer(CMonomer* monomer) { _monomers.emplace_back(unique_ptr<CMonomer>(monomer));}
    /// Get monomer at an index
    /// @note no check on index
    CMonomer* getCMonomer(int index) {return _monomers[index].get();}
    
    ///Get list of reactions associated
    const unordered_set<ReactionBase*>& getInternalReactions() {return _internalReactions;}
    const unordered_set<CCylinder*>& getReactingCylinders() {return _reactingCylinders;}
    ///Get map of reactions associated
    unordered_map<CCylinder*, unordered_set<ReactionBase*>>& getCrossCylinderReactions() {return _crossCylinderReactions;}
    
    //@{
    /// Reaction management function
    void addInternalReaction(ReactionBase* r);
    void addCrossCylinderReaction(CCylinder* other, ReactionBase* r);
    void addReactingCylinder(CCylinder* other);

    void removeInternalReaction(ReactionBase* r);
    void removeAllInternalReactions();
    
    void removeCrossCylinderReactions(CCylinder* other);
    void removeAllCrossCylinderReactions();
    
    void removeReactingCylinder(CCylinder* other);
    void removeAllReactingCylinders();
    //@}
    
    /// Print
    virtual void printCCylinder();
    
    /// Get size in number of monomers
    short getSize() {return _size;}
    
};


#endif /* defined(__CytoSim__CCylinder__) */
