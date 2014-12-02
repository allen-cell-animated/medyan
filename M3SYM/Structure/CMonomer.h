
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

#ifndef M3SYM_CMonomer_h
#define M3SYM_CMonomer_h

#include "common.h"

#include "Species.h"
#include "Compartment.h"

/// Rpresents a container for all [Species] (@ref Species) that could be contained in a
/// particular filament element at a given position.
/*!
 *  CMonomer provides a container to hold all [Species] (@ref Species) that are possibly held at a given
 *  filament position. The species are held in an standard vector. Functions to lookup species
 *  as well as a filament element checker are provided.
 */
class CMonomer {
    
    //@{
    /// Species array
    SpeciesFilament** _speciesFilament;
    SpeciesPlusEnd**  _speciesPlusEnd;
    SpeciesMinusEnd** _speciesMinusEnd;
    SpeciesBound**    _speciesBound;
    SpeciesLinker**   _speciesLinker;
    SpeciesMotor**    _speciesMotor;
    //@}
    
public:
    /// Constructor does nothing but memset arrays
    CMonomer();
    
    /// Default destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~CMonomer () noexcept;
    
    /// Copy constructor
    /// This constructor will create a new CMonomer, identical to the copied, in a new compartment. The
    /// original species will remain intact, while new identical species will be initialized.
    CMonomer(const CMonomer& rhs, Compartment* c);
    /// Assignment is not allowed
    CMonomer& operator=(CMonomer &rhs) = delete;
    
    /// Clone a CMonomer. Transfers all [CBounds] (@ref CBound) to the new CMonomer.
    virtual CMonomer* clone(Compartment* c) {
        return new CMonomer(*this, c);
    }
    
    ///@{
    /// Add a species
    /// @note should only be called at initialization
    void addSpeciesFilament (SpeciesFilament* s);
    void addSpeciesPlusEnd  (SpeciesPlusEnd* s);
    void addSpeciesMinusEnd (SpeciesMinusEnd* s);
    void addSpeciesBound    (SpeciesBound* s);
    void addSpeciesLinker   (SpeciesLinker* s);
    void addSpeciesMotor    (SpeciesMotor* s);
    //@}
    
    ///Print the species in this filament element
    void print();

    //@{
    /// Get species at a specific index
    /// @note no check on this index. The index value of a species is stored in the chemical initializer
    /// when all reactions are initialized from the chemical input file
    SpeciesFilament* speciesFilament (int index)  {return _speciesFilament[index];}
    SpeciesPlusEnd*  speciesPlusEnd  (int index)  {return _speciesPlusEnd[index];}
    SpeciesMinusEnd* speciesMinusEnd (int index)  {return _speciesMinusEnd[index];}
    SpeciesBound*    speciesBound    (int index)  {return _speciesBound[index];}
    SpeciesLinker*   speciesLinker   (int index)  {return _speciesLinker[index];}
    SpeciesMotor*    speciesMotor    (int index)  {return _speciesMotor[index];}
    //@}
    
    ///Check if this filament element is valid. Involves checking copy numbers
    virtual bool checkSpecies(int sum) {return true;}
    
};




#endif /* defined(__CytoSim__CMonomer__) */
