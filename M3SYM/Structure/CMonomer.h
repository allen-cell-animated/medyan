
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

#ifndef M3SYM_CMonomer_h
#define M3SYM_CMonomer_h

#include "common.h"

#include "Species.h"
#include "Compartment.h"

/// Represents a container for all Species that could be contained in a
/// particular filament element at a given position.
/*!
 *  CMonomer provides a container to hold all Species that are possibly held at a given
 *  filament position. The species are held in an standard vector.
 */
class CMonomer {
    
    //@{
    /// Species array
    SpeciesFilament** _speciesFilament;
    SpeciesPlusEnd**  _speciesPlusEnd;
    SpeciesMinusEnd** _speciesMinusEnd;
    
    SpeciesBound**  _speciesBound;
    SpeciesLinker** _speciesLinker;
    SpeciesMotor**  _speciesMotor;
    SpeciesBrancher** _speciesBrancher;
    //@}
    
public:
    /// Constructor does nothing but memset arrays
    CMonomer();
    
    /// Default destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is
    /// a gcc bug (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~CMonomer () noexcept;
    
    /// Copy constructor
    /// This constructor will create a new CMonomer, identical to the copied, in a new
    /// compartment. The original species will remain intact, while new identical
    /// species will be initialized.
    CMonomer(const CMonomer& rhs, Compartment* c);
    /// Assignment is not allowed
    CMonomer& operator=(CMonomer &rhs) = delete;
    
    /// Clone a CMonomer. Transfers all [CBounds] (@ref CBound) to the new CMonomer.
    virtual CMonomer* clone(Compartment* c) {
        return new CMonomer(*this, c);
    }
    
    ///@{
    /// Add Species
    /// @note should only be called at initialization
    void addSpeciesFilament (SpeciesFilament* s);
    void addSpeciesPlusEnd  (SpeciesPlusEnd* s);
    void addSpeciesMinusEnd (SpeciesMinusEnd* s);
    void addSpeciesBound    (SpeciesBound* s);
    void addSpeciesLinker   (SpeciesLinker* s);
    void addSpeciesMotor    (SpeciesMotor* s);
    void addSpeciesBrancher (SpeciesBrancher* s);
    //@}
    
    ///Print the Species
    void print();

    //@{
    /// Get Species at a specific index
    /// @note no check on this index. The index value of a species is stored in the
    /// chemical initializer when all reactions are initialized from the chemical input
    /// file.
    inline SpeciesFilament* speciesFilament(int index) {return _speciesFilament[index];}
    inline SpeciesPlusEnd*  speciesPlusEnd (int index) {return _speciesPlusEnd[index];}
    inline SpeciesMinusEnd* speciesMinusEnd(int index) {return _speciesMinusEnd[index];}
    
    inline SpeciesBound*    speciesBound   (int index) {return _speciesBound[index];}
    inline SpeciesLinker*   speciesLinker  (int index) {return _speciesLinker[index];}
    inline SpeciesMotor*    speciesMotor   (int index) {return _speciesMotor[index];}
    inline SpeciesBrancher* speciesBrancher(int index) {return _speciesBrancher[index];}
    //@}
    
    //@{
    /// Get the active species of the corresponding type
    /// @note these functions will return -1 if the corresponding type has no marked
    /// species at this point.
    short activeSpeciesFilament();
    short activeSpeciesPlusEnd();
    short activeSpeciesMinusEnd();
    
    short activeSpeciesBound();
    short activeSpeciesLinker();
    short activeSpeciesMotor();
    short activeSpeciesBrancher();
    //@
};

#endif
