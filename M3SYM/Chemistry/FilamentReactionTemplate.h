
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

#ifndef M3SYM_FilamentReactionTemplate_h
#define M3SYM_FilamentReactionTemplate_h

#include <vector>
#include <cmath>

#include "common.h"

#include "Species.h"

#include "SysParams.h"

///Enumeration for direction of reaction
enum FilamentReactionDirection {
    FORWARD, BACKWARD, INPLACE
};

//FORWARD DECLARATIONS
class SubSystem;
class CCylinder;

/// To store Filament chemical reaction information read from an input file
/*!
 *  FilamentReactionTemplate is used to store a filament reaction. It contains vectors
 *  of tuples that represent the position in the CMonomer in which the species is stored 
 *  (for products and reactants), as well as the rate of the reaction. The integer value 
 *  that is the position of the species in the CMonomer vector is held by the 
 *  ChemManager.
 *  @note if the species is a bulk or diffusing species, the integer molecule value in 
 *  the tuple stored in the SpeciesNamesDB.
 *
 *  This class also has functions to add the filament reaction to a CCylinder, as well 
 *  as add a connection reaction between two neighboring [CCylinders](@ref CCylinder).
 */
class FilamentReactionTemplate {
    
    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to initialize and
                           ///< call chemical callbacks

    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    /// Info about binding chemistry
    short _empty = 0;
    
    float _rate; ///< Rate of reaction
    
public:
    FilamentReactionTemplate(vector<tuple<int, SpeciesType>> reactants,
                             vector<tuple<int, SpeciesType>> products,
                             float rate)
        : _reactants(reactants), _products(products), _rate(rate) {

#if !defined(REACTION_SIGNALING)
        cout << "Any filament reaction relies on reaction signaling. Please"
            << " set this compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
    }
    ~FilamentReactionTemplate() {}

    /// Add this chemical reaction. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc) = 0;
    
    /// Add this chemical reaction along a filament
    /// @note assumes cc1 and cc2 are in order, that is, cc2 is the next
    /// cylinder after cc1
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
};

/// Manager for polymerization at plus end of Filament
class PolyPlusEndManager : public FilamentReactionTemplate {
    
public:
    PolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~PolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for polymerization at minus end of Filament
class PolyMinusEndManager : public FilamentReactionTemplate {
    
public:
    PolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~PolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};


/// Manager for depolymerization at plus end of Filament
class DepolyPlusEndManager : public FilamentReactionTemplate {
    
public:
    DepolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~DepolyPlusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for depolymerization at minus end of Filament
class DepolyMinusEndManager : public FilamentReactionTemplate {
    
public:
    DepolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                          vector<tuple<int, SpeciesType>> products,
                          float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~DepolyMinusEndManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost walking
class MotorWalkFManager : public FilamentReactionTemplate {
    
public:
    ///default constructor and destructor
    MotorWalkFManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~MotorWalkFManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Manager for MotorGhost walking
class MotorWalkBManager : public FilamentReactionTemplate {
    
public:
    ///default constructor and destructor
    MotorWalkBManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~MotorWalkBManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


/// Manager for Filament aging
class AgingManager : public FilamentReactionTemplate {
    
public:
    AgingManager(vector<tuple<int, SpeciesType>> reactants,
                 vector<tuple<int, SpeciesType>> products,
                 float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~AgingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Manager for Filament destruction
class DestructionManager : public FilamentReactionTemplate {
    
    
public:
    DestructionManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~DestructionManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
    
};

/// Manager for severing a Filament
class SeveringManager : public FilamentReactionTemplate {
    
public:
    SeveringManager(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~SeveringManager() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};


//@{
/// Helper for tuple getter
inline int getInt(tuple<int, SpeciesType> x) {return get<0>(x);}
inline SpeciesType getType(tuple<int, SpeciesType> x) {return get<1>(x);}
//@}

#endif
