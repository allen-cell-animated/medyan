
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

#ifndef M3SYM_ReactionManager_h
#define M3SYM_ReactionManager_h

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
 *  FilamentRxnManager is used to store a filament reaction. It contains vectors
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
class FilamentRxnManager {
    
    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to initialize and
                           ///< call chemical callbacks

    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    /// Info about binding chemistry
    short _empty = 0;
    
    float _rate; ///< Rate of reaction
    
    vector<short> _bindingSites; ///< The binding sites on the Cylinder
    
public:
    FilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate)
        : _reactants(reactants), _products(products), _rate(rate) {

#if !defined(REACTION_SIGNALING)
        cout << "Any filament reaction relies on reaction signaling. Please"
            << " set this compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
        //Figure out the binding sites
        int deltaBinding = SysParams::Geometry().cylinderIntSize /
        SysParams::Chemistry().numBindingSites;
        
        int firstBindingSite = deltaBinding / 2 + 1;
        int bindingCount = firstBindingSite;
        
        //add all other binding sites
        while(bindingCount < SysParams::Geometry().cylinderIntSize) {
            
            _bindingSites.push_back(bindingCount);
            bindingCount += deltaBinding;
        }
    }
    ~FilamentRxnManager() {}

    /// Add this chemical reaction. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc) = 0;
    
    /// Add this chemical reaction along a filament
    /// @note assumes cc1 and cc2 are in order, that is, cc2 is the next
    /// cylinder after cc1
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
};

/// To store binding reactions, including Linker, MotorGhost, and BranchingPoint binding

/*!
 *  BindingRxnManager is used to store a binding reaction on filaments in compartments.
 *  Contains the binding reaction, possible binding sites, and integers representing the 
 *  binding species involved in the reaction.
 *  Classes that extend this will implement their own data structures for holding 
 *  possible reaction sites, etc.
 *  The main function of this class is to call updatePossibleBindings(), which will 
 *  update the possible binding sites if the binding reaction is called in this compartment.
 */
class BindingRxnManager {

    friend class SimpleManagerImpl;
    
protected:
    ReactionBase* _bindingReaction; ///< The binding reaction for this compartment
    
    vector<short> _bindingSites; ///< The binding sites on the Cylinder
    
    //@{
    /// Info about binding chemistry
    short _empty = 0;
    short _bound;
    //@}

public:
    BindingRxnManager(ReactionBase* reaction, short bound)
    
        : _bindingReaction(reaction), _bound(bound) {
            
#if !defined(REACTION_SIGNALING)
            cout << "Any filament binding reaction relies on reaction signaling. Please"
            << " set this compilation macro and try again. Exiting." << endl;
            exit(EXIT_FAILURE);
#endif
        //Figure out the binding sites
        int deltaBinding = SysParams::Geometry().cylinderIntSize /
                           SysParams::Chemistry().numBindingSites;
                                  
        int firstBindingSite = deltaBinding / 2 + 1;
        int bindingCount = firstBindingSite;
        
        //add all other binding sites
        while(bindingCount < SysParams::Geometry().cylinderIntSize) {
            
            _bindingSites.push_back(bindingCount);
            bindingCount += deltaBinding;
        }
    }
    
    ///update the possible binding reactions that could occur
    virtual void updatePossibleBindings(CCylinder* cc) = 0;
};

//@{
/// Helper for tuple getter
inline int getInt(tuple<int, SpeciesType> x) {return get<0>(x);}
inline SpeciesType getType(tuple<int, SpeciesType> x) {return get<1>(x);}
//@}

#endif
