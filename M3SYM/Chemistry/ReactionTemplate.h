
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

#ifndef M3SYM_ReactionTemplate_h
#define M3SYM_ReactionTemplate_h

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
 *  ChemTemplate.
 *  @note if the species is a bulk or diffusing species, the integer molecule value in 
 *  the tuple stored in the SpeciesNamesDB.
 *
 *  This class also has functions to add the filament reaction to a CCylinder, as well 
 *  as add a connection reaction between two neighboring [CCylinders](@ref CCylinder).
 */
class FilamentReactionTemplate {
    
friend class CController;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to initialize and
                           ///< call chemical callbacks

    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    float _rate; ///< Rate of reaction
    
public:
    FilamentReactionTemplate(vector<tuple<int, SpeciesType>> reactants,
                             vector<tuple<int, SpeciesType>> products,
                             float rate)
        : _reactants(reactants), _products(products), _rate(rate) {

#if !defined(REACTION_SIGNALING) || !defined(RSPECIES_SIGNALING)
        cout << "Any filament reaction relies on reaction and species signaling. Please"
            << " set this compilation macros and try again. Exiting." << endl;
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

/// Template for polymerization at plus end of Filament
class PolyPlusEndTemplate : public FilamentReactionTemplate {
    
public:
    PolyPlusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~PolyPlusEndTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Template for polymerization at minus end of Filament
class PolyMinusEndTemplate : public FilamentReactionTemplate {
    
public:
    PolyMinusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~PolyMinusEndTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};


/// Template for depolymerization at plus end of Filament
class DepolyPlusEndTemplate : public FilamentReactionTemplate {
    
public:
    DepolyPlusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~DepolyPlusEndTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Template for depolymerization at minus end of Filament
class DepolyMinusEndTemplate : public FilamentReactionTemplate {
    
public:
    DepolyMinusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                          vector<tuple<int, SpeciesType>> products,
                          float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~DepolyMinusEndTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Template for MotorGhost walking
class MotorWalkFTemplate : public FilamentReactionTemplate {
    
public:
    ///default constructor and destructor
    MotorWalkFTemplate(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~MotorWalkFTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

/// Template for MotorGhost walking
class MotorWalkBTemplate : public FilamentReactionTemplate {
    
public:
    ///default constructor and destructor
    MotorWalkBTemplate(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~MotorWalkBTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


/// Template for Filament aging
class AgingTemplate : public FilamentReactionTemplate {
    
public:
    AgingTemplate(vector<tuple<int, SpeciesType>> reactants,
                 vector<tuple<int, SpeciesType>> products,
                 float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~AgingTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};

/// Template for Filament destruction
class DestructionTemplate : public FilamentReactionTemplate {
    
    
public:
    DestructionTemplate(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~DestructionTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
    
};

/// Template for severing a Filament
class SeveringTemplate : public FilamentReactionTemplate {
    
public:
    SeveringTemplate(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate)
    : FilamentReactionTemplate(reactants, products, rate) {}
    ~SeveringTemplate() {}
    
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {}
};


//@{
/// Helper for tuple getter
inline int getInt(tuple<int, SpeciesType> x) {return get<0>(x);}
inline SpeciesType getType(tuple<int, SpeciesType> x) {return get<1>(x);}
//@}

#endif
