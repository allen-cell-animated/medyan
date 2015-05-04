
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_ReactionManagerImpl_h
#define M3SYM_ReactionManagerImpl_h

#include <unordered_set>

#include "common.h"

#include "ReactionManager.h"
#include "NeighborListContainer.h"

/// Manager for polymerization at plus end of Filament
class PolyPlusEndManager : public FilamentRxnManager {
    
public:
    PolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~PolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for polymerization at minus end of Filament
class PolyMinusEndManager : public FilamentRxnManager {
    
public:
    PolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~PolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};


/// Manager for depolymerization at plus end of Filament
class DepolyPlusEndManager : public FilamentRxnManager {
    
public:
    DepolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~DepolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for depolymerization at minus end of Filament
class DepolyMinusEndManager : public FilamentRxnManager {
    
public:
    DepolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                          vector<tuple<int, SpeciesType>> products,
                          float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~DepolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost walking
class MotorWalkFManager : public FilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkFManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkFManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost walking
class MotorWalkBManager : public FilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkBManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkBManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


/// Manager for Filament aging
class AgingManager : public FilamentRxnManager {
    
public:
    AgingManager(vector<tuple<int, SpeciesType>> reactants,
                 vector<tuple<int, SpeciesType>> products,
                 float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~AgingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for Filament destruction
class DestructionManager : public FilamentRxnManager {
    
    
public:
    DestructionManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~DestructionManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
    
};

/// Manager for severing a Filament
class SeveringManager : public FilamentRxnManager {
    
public:
    SeveringManager(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate)
    : FilamentRxnManager(reactants, products, rate) {}
    ~SeveringManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for Filament and BranchPoint creation
class BranchingManager : public BindingRxnManager {
    
private:
    ///possible bindings at current state
    unordered_set<tuple<CCylinder*, short>> _possibleBindings;
    
public:
    BranchingManager(ReactionBase* reaction, short bound)
        : BindingRxnManager(reaction, bound) {}
    
    ~BranchingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc);
};

/// Manager for Linker binding.
/*!
 *  LinkerBindingManager controls linker binding in a compartment.
 *  Also a subclass of CCNLContainer, contains a cylinder neighbors list of
 *  cylinders within range of binding for this reaction
 */
class LinkerBindingManager : public BindingRxnManager, public CCNLContainer {
    
private:
    
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_map<tuple<CCylinder*, short>, tuple<Cylinder*, short>> _possibleBindings;
    
public:
    LinkerBindingManager(ReactionBase* reaction, short bound, float rMax, float rMin)
    
    : BindingRxnManager(reaction, bound),
      CCNLContainer(rMax + SysParams::Geometry().cylinderSize,
                max(rMin - SysParams::Geometry().cylinderSize, 0.0)),
    
     _rMin(rMin), _rMax(_rMax) {}
    
    ~LinkerBindingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc);
};

/// Manager for MotorGhost binding
/*!
 *  MotorBindingManager controls motor binding in a compartment.
 *  Also a subclass of CCNLContainer, contains a cylinder neighbors list of
 *  cylinders within range of binding for this reaction
 */
class MotorBindingManager : public BindingRxnManager, public CCNLContainer {
    
private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_map<tuple<CCylinder*, short>, tuple<Cylinder*, short>> _possibleBindings;
    
public:
    MotorBindingManager(ReactionBase* reaction, short bound, float rMax, float rMin)
    
    : BindingRxnManager(reaction, bound),
      CCNLContainer(rMax + SysParams::Geometry().cylinderSize,
                max(rMin - SysParams::Geometry().cylinderSize, 0.0)),
    
    _rMin(rMin), _rMax(_rMax) {}
    
    ~MotorBindingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc);
};


#endif
