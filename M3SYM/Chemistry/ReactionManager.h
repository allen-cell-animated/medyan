
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

#ifndef M3SYM_ReactionTemplate_h
#define M3SYM_ReactionTemplate_h

#include <vector>

#include "common.h"

#include "NeighborListContainer.h"
#include "CompartmentContainer.h"

///Enumeration for direction of reaction
enum FilamentReactionDirection {
    FORWARD, BACKWARD, INPLACE
};

//FORWARD DECLARATIONS
class SubSystem;
class CCylinder;

///InternalFilamentRxnManager is a class to store filament chemical reaction information read from an input file
/*!
 *  InternalFilamentRxnManager is used to store a filament reaction. It contains vectors of tuples that represent
 *  the position in the [CMonomer] (@ref CMonomer) in which the species is stored (for products and reactants), as well as the rate
 *  of the reaction. The integer value that is the position of the species in the CMonomer vector is held by the [ChemManager] (@ref ChemManager).
 *
 *  @note if the species is a bulk or diffusing species, the integer molecule value in the tuple 
 *  stored in the [SpeciesNamesDB] (@ref SpeciesNamesDB).
 *
 *  This class also has functions to add the filament reaction to a [CCylinder] (@ref CCylinder), as well as add a 
 *  connection reaction between two neighboring CCylinders.
 */
class InternalFilamentRxnManager {
    
    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to initialize and call chemical callbacks
    
    ///Species identifier vectors
    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    float _rate; ///< Rate of reaction
    
public:
    InternalFilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                               vector<tuple<int, SpeciesType>> products,
                               float rate) : _reactants(reactants), _products(products), _rate(rate) {}
    ~InternalFilamentRxnManager() {}

    ///Add this chemical reaction to a CCylinder. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc) = 0;
    
    ///Add this chemical reaction cross two CCylinders.
    ///@note assumes cc1 and cc2 are in order, that is, cc2 is the next cylinder after cc1 
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
};

/// Manager for polymerization at plus end
class PolyPlusEndManager : public InternalFilamentRxnManager {

public:
    PolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~PolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for polymerization at minus end
class PolyMinusEndManager : public InternalFilamentRxnManager {
    
public:
    PolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~PolyMinusEndManager() {}

    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


/// Manager for depolymerization at plus end
class DepolyPlusEndManager : public InternalFilamentRxnManager {
    
public:
    DepolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DepolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for depolymerization at minus end
class DepolyMinusEndManager : public InternalFilamentRxnManager {
    
public:
    DepolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                          vector<tuple<int, SpeciesType>> products,
                          float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DepolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for basic binding
class BasicBindingManager : public InternalFilamentRxnManager {
    
public:
    BasicBindingManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~BasicBindingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {};
};

/// Manager for all unbinding
class UnbindingManager : public InternalFilamentRxnManager {
    
public:
    UnbindingManager(vector<tuple<int, SpeciesType>> reactants,
                     vector<tuple<int, SpeciesType>> products,
                     float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~UnbindingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {};

};

/// Manager for motor walking
class MotorWalkFManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkFManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkFManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);

};

/// Manager for motor walking
class MotorWalkBManager : public InternalFilamentRxnManager {
    
public:
    MotorWalkBManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkBManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
    
};

/// CrossFilamentRxnManager is a class to store cross-filament reactions, including linker and motor binding

/*!
 *  CrossFilamentRxnManager is used to store a cross-filament reaction. It contains vectors of tuples that represent
 *  the position in the [CMonomer] (@ref CMonomer) in which the species is stored (for products and reactants), as well as the rate
 *  of the reaction and direction. The integer value that is the position of the species in the CMonomer vector
 *  is held by the [ChemManager] (@ref ChemManager). Also contains the range of this reaction.
 *
 *  Also a subclass of [CylinderNLContainer] (@ref CylinderNLContainer), contains a cylinder neighbors list of active reactions
 *
 *  @note if the species is a bulk or diffusing species, the integer molecule value in the tuple stored in the [SpeciesNamesDB] (@ref SpeciesNamesDB).
 *
 *  This class also has functions to add the cross filament reaction to two [CCylinders] (@ref CCylinder)
 */

class CrossFilamentRxnManager : public CylinderNLContainer {

    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to intialize and call chemical callbacks
    
    ///Species identifier vectors
    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    float _rate; ///< Rate of reaction
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    int _reactionID; ///Unique ID of this reaction

public:
    CrossFilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                            vector<tuple<int, SpeciesType>> products,
                            float rate, float rMin, float rMax, int ID)
                            : CylinderNLContainer(rMax, rMin, true),
                              _reactants(reactants), _products(products), _rate(rate), _rMin(rMin), _rMax(rMax), _reactionID(ID) {}
    ~CrossFilamentRxnManager() {}
    
    /// Add this chemical reaction to two ccylinders if within the reaction range
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
    
    //@{
    /// Getter for parameters of reaction
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    int getReactionID() {return _reactionID;}
    //@}
    
};

/// Manager for linker binding
class LinkerBindingManager : public CrossFilamentRxnManager {
    
public:
    LinkerBindingManager(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate, float rMin, float rMax, int ID)
                    : CrossFilamentRxnManager(reactants, products, rate, rMin, rMax, ID) {}
    ~LinkerBindingManager() {}

    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for motor binding
class MotorBindingManager : public CrossFilamentRxnManager {
    
public:
    MotorBindingManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate, float rMin, float rMax, int ID)
                         : CrossFilamentRxnManager(reactants, products, rate, rMin, rMax, ID) {}
    ~MotorBindingManager() {}
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};



#endif
