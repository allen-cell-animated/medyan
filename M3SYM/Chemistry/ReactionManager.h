
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

#include "NeighborListContainer.h"
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
 *  InternalFilamentRxnManager is used to store a filament reaction. It contains vectors 
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
class InternalFilamentRxnManager {
    
    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to initialize and
                           ///< call chemical callbacks

    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    float _rate; ///< Rate of reaction
    
    vector<short> _bindingSites; ///< Binding sites on Cylinder
    
public:
    InternalFilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                               vector<tuple<int, SpeciesType>> products,
                               float rate)
        : _reactants(reactants), _products(products), _rate(rate) {
    
        
#if !defined(REACTION_SIGNALING)
        cout << "Any filament-related reaction relies on reaction signaling. Please"
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
    ~InternalFilamentRxnManager() {}

    /// Add this chemical reaction. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc) = 0;
    
    /// Add this chemical reaction along a filament
    /// @note assumes cc1 and cc2 are in order, that is, cc2 is the next
    /// cylinder after cc1
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
};

/// To store cross-filament reactions, including Linker and MotorGhost binding

/*!
 *  CrossFilamentRxnManager is used to store a cross-filament reaction. It contains 
 *  vectors of tuples that represent the position in the CMonomer in which the species 
 *  is stored (for products and reactants), as well as the rate of the reaction and 
 *  direction. The integer value that is the position of the species in the CMonomer 
 *  vector is held by the ChemManager. Also contains the range of this reaction.
 *  Also a subclass of CCNLContainer, contains a cylinder neighbors list of
 *  active reactions.
 *  @note if the species is a bulk or diffusing species, the integer molecule value in 
 *  the tuple stored in the SpeciesNamesDB.
 *  This class also has functions to add the cross filament reaction to two 
 *  [CCylinders] (@ref CCylinder).
 */
class CrossFilamentRxnManager : public CCNLContainer {

    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps; ///< A subsystem pointer to intialize
                           ///< and call chemical callbacks
    
    ///Species identifier vectors
    vector<tuple<int,SpeciesType>> _reactants; ///< Reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< Products in this reaction
    
    float _onRate; ///< On-rate for interaction
    float _offRate; ///< Off-rate for interaction
    
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    vector<short> _bindingSites; ///< The binding sites on the Cylinder

public:
    CrossFilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                            vector<tuple<int, SpeciesType>> products,
                            float onRate, float offRate, float rMax, float rMin)
    
        : CCNLContainer(rMax + SysParams::Geometry().cylinderSize,
                        max(rMin - SysParams::Geometry().cylinderSize, 0.0), true),
    
        _reactants(reactants), _products(products),
        _onRate(onRate), _offRate(offRate), _rMin(rMin), _rMax(rMax) {
            
#if !defined(REACTION_SIGNALING)
            cout << "Any filament-related reaction relies on reaction signaling. Please"
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
    ~CrossFilamentRxnManager() {}
    
    /// Add this chemical reaction to two ccylinders if within the reaction range
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
    
    //@{
    /// Getter for parameters of reaction
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    //@}
    
};

//@{
/// Helper for tuple getter
inline int getInt(tuple<int, SpeciesType> x) {return get<0>(x);}
inline SpeciesType getType(tuple<int, SpeciesType> x) {return get<1>(x);}
//@}

#endif
