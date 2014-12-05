
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

#ifndef M3SYM_CBound_h
#define M3SYM_CBound_h

#include "common.h"

#include "Species.h"

//FORWARD DECLARATIONS
class Compartment;

/// Represents a chemical object that is bound to a Filament.
/*!
 *  The CBound class is an abstract representation of a chemically bound object to a Filament
 *  (Could be a Linker, MotorGhost, BranchingPoint, etc). Each CBound object has a pointer
 *  to the corresponding [SpeciesBound] (@ref SpeciesBound) on a Filament. Different implementations of CBound class will
 *  have different functions to bind, move, etc. See documentation of subclass for more details on function.
 */
class CBound {
    
protected:
    SpeciesBound* _firstSpecies = nullptr; ///< Corresponding first species on Filament
    SpeciesBound* _secondSpecies = nullptr; ///< Corresponding second species on Filament
    
    Compartment* _compartment; ///< Compartment this CBound is in
    
public:
    /// Constructor, just sets species
    CBound(Compartment* c) : _compartment(c) {}
    
    /// Virtual destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~CBound() noexcept {
        _firstSpecies->removeCBound();
        _secondSpecies->removeCBound();
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
    
    /// Get compartment that this CBound is in
    Compartment* getCompartment() {return _compartment;}
};

#endif
