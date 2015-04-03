
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

#ifndef M3SYM_CBound_h
#define M3SYM_CBound_h

#include "common.h"

#include "Species.h"
#include "ReactionBase.h"

//FORWARD DECLARATIONS
class Compartment;
class SubSystem;
class CCylinder;

/// Represents a chemical object that is bound to a Filament.
/*!
 *  The CBound class is an abstract representation of a chemically bound object to a 
 *  Filament (Could be a Linker, MotorGhost, BranchingPoint, etc). Each CBound object 
 *  has a pointer to the corresponding [SpeciesBound] (@ref SpeciesBound) on a Filament. 
 *  Different implementations of CBound class will have different functions to bind,
 *  move, etc. See documentation of subclass for more details on function.
 */
class CBound {
    
protected:
    SpeciesBound* _firstSpecies = nullptr; ///< Corresponding first species on Filament
    SpeciesBound* _secondSpecies = nullptr;///< Corresponding second species on Filament
    
    Compartment* _compartment; ///< Compartment this CBound is in
    
    CCylinder* _cc1 = nullptr; ///< Pointer to first CCylinder
    CCylinder* _cc2 = nullptr; ///< Pointer to second CCylinder
    
    //@{
    ///Reaction rates
    float _onRate = 0.0;
    float _offRate = 0.0;
    //@}
    
    ReactionBase* _offRxn; ///< The off reaction for this bound object
    
public:
    /// Constructor, just sets species
    CBound(Compartment* c, CCylinder* cc1, CCylinder* cc2)
        : _compartment(c), _cc1(cc1), _cc2(cc2) {}
    
    /// Virtual destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~CBound() noexcept {
        if(_firstSpecies != nullptr) _firstSpecies->removeCBound();
        if(_secondSpecies != nullptr) _secondSpecies->removeCBound();
    }
    
    /// Set first species
    void setFirstSpecies(SpeciesBound* species) {
        ///remove from old first species
        if(_firstSpecies != nullptr) _firstSpecies->removeCBound();
        
        _firstSpecies = species;
        _firstSpecies->setCBound(this);
    }
    /// Get first species
    SpeciesBound* getFirstSpecies() {return _firstSpecies;}
    
    /// Set second species
    void setSecondSpecies(SpeciesBound* species) {
        ///remove from old second species
        if(_secondSpecies != nullptr) _secondSpecies->removeCBound();
        
        _secondSpecies = species;
        _secondSpecies->setCBound(this);
    }
    /// Get second species
    SpeciesBound* getSecondSpecies() {return _secondSpecies;}
    
    /// Set first CCylinder
    void setFirstCCylinder(CCylinder* cc) {_cc1 = cc;}
    /// Get first CCylinder
    CCylinder* getFirstCCylinder() {return _cc1;}
    
    /// Set second CCylinder
    void setSecondCCylinder(CCylinder* cc) {_cc2 = cc;}
    /// Get first CCylinder
    CCylinder* getSecondCCylinder() {return _cc2;}
    
    /// Get compartment that this CBound is in
    Compartment* getCompartment() {return _compartment;}
    
    //@{
    /// On rate management
    void setOnRate(float rate) {_onRate = rate;}
    double getOnRate(){return _onRate;}
    //@}
    
    //@{
    /// Off rate management
    void setOffRate(float rate) {_offRate = rate;}
    double getOffRate(){return _offRate;}
    //@}
    
    
    //@{
    /// Off reaction management
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps) = 0;
    
    void setOffReaction(ReactionBase* offRxn) {
        _offRxn = offRxn;
        _offRxn->setCBound(this);
    }
    ReactionBase* getOffReaction() {return _offRxn;}
    //@}
};

#endif
