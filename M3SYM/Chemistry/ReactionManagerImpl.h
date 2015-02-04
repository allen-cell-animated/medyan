
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

#include "common.h"

#include "ReactionManager.h"

/// Manager for polymerization at plus end of Filament
class PolyPlusEndManager : public InternalFilamentRxnManager {
    
public:
    PolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~PolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for polymerization at minus end of Filament
class PolyMinusEndManager : public InternalFilamentRxnManager {
    
public:
    PolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~PolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};


/// Manager for depolymerization at plus end of Filament
class DepolyPlusEndManager : public InternalFilamentRxnManager {
    
public:
    DepolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DepolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for depolymerization at minus end of Filament
class DepolyMinusEndManager : public InternalFilamentRxnManager {
    
public:
    DepolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                          vector<tuple<int, SpeciesType>> products,
                          float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DepolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost walking
class MotorWalkFManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkFManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkFManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost walking
class MotorWalkBManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkBManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkBManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


/// Manager for Filament aging
class AgingManager : public InternalFilamentRxnManager {
    
public:
    AgingManager(vector<tuple<int, SpeciesType>> reactants,
                 vector<tuple<int, SpeciesType>> products,
                 float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~AgingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for Filament destruction
class DestructionManager : public InternalFilamentRxnManager {
    
    
public:
    DestructionManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DestructionManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
    
};

/// Manager for severing a Filament
class SeveringManager : public InternalFilamentRxnManager {
    
public:
    SeveringManager(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate)
    : InternalFilamentRxnManager(reactants, products, rate) {}
    ~SeveringManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for Filament and BranchPoint creation
class BranchingManager : public InternalFilamentRxnManager {
    
private:
    float _offRate; ///< Off rate for this branchingpoint
    
public:
    BranchingManager(vector<tuple<int, SpeciesType>> reactants,
                     vector<tuple<int, SpeciesType>> products,
                     float rate, float offRate)
    : InternalFilamentRxnManager(reactants, products, rate), _offRate(offRate) {}
    ~BranchingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for Linker binding and unbinding
class LinkerRxnManager : public CrossFilamentRxnManager {
    
public:
    LinkerRxnManager(vector<tuple<int, SpeciesType>> reactants,
                     vector<tuple<int, SpeciesType>> products,
                     float onRate, float offRate, float rMax, float rMin)
    : CrossFilamentRxnManager(reactants, products, onRate, offRate, rMin, rMax) {}
    ~LinkerRxnManager() {}
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost binding and unbinding
class MotorRxnManager : public CrossFilamentRxnManager {
    
public:
    MotorRxnManager(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float onRate, float offRate, float rMax, float rMin)
    : CrossFilamentRxnManager(reactants, products, onRate, offRate, rMin, rMax) {}
    ~MotorRxnManager() {}
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


#endif
