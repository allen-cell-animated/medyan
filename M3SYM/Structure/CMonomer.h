
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
 *  CMonomer provides a container to hold all Species that are held at a given
 *  filament position. The species are held in an standard vector.
 */
class CMonomer {
    
friend class ChemManager;
friend class CCylinder;

private:
    //@{
    /// Species array
    SpeciesFilament** _speciesFilament;
    SpeciesBound** _speciesBound;
    //@}
    
    //@{
    /// Species index vectors
    /// These index vectors are used to access the correct species in the actual
    /// species arrays. Each index in the species index vector corresponds to an
    /// offset for that species in the species array.
    static vector<short> _speciesFilamentIndex;
    static vector<short> _speciesBoundIndex;
    //@}
    
    ///Number of species for each type
    static short _numFSpecies;
    static short _numBSpecies;
    
public:
    /// Constructor does nothing but memset arrays
    CMonomer();
    
    /// Default destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL containers.
    /// This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly be fixed in the future.
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
    
    ///Print the Species
    void print();

    //@{
    /// Get Species at a specific index
    /// @note no check on this index. The index value of a species is stored in the
    /// chemical initializer when all reactions are initialized from the chemical input
    /// file.
    SpeciesFilament* speciesFilament (int index);
    SpeciesFilament* speciesPlusEnd  (int index);
    SpeciesFilament* speciesMinusEnd (int index);
    
    SpeciesBound* speciesBound    (int index);
    SpeciesBound* speciesLinker   (int index);
    SpeciesBound* speciesMotor    (int index);
    SpeciesBound* speciesBrancher (int index);
    //@}
    
    //@{
    /// Get the active species of the corresponding type
    /// @note these functions will return -1 if the corresponding type has no marked
    /// species at this point.
    short activeSpeciesFilament();
    short activeSpeciesPlusEnd();
    short activeSpeciesMinusEnd();
    
    short activeSpeciesLinker();
    short activeSpeciesMotor();
    short activeSpeciesBrancher();
    //@
    
    /// Check the consistency of the CMonomer for debugging.
    bool isConsistent();
    
};

#endif
