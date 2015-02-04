
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

#ifndef M3SYM_MLinker_h
#define M3SYM_MLinker_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents the mechanical component of a Linker.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a Linker,
 *  and its associated equilibrium constants. Initial length of a Linker is determined
 *  by the condition of zero initial stress, i.e., it calculated within the constructor
 *  at initiation. A Linker heads positions on a segment (between two consecutive beads 
 *  on a Filament) determined by two numbers (double from 0 to 1) position1 and 
 *  position2.
 */
class MLinker {
    
public:
    double stretchForce = 0.0; ///< Stretching force of linker at current state
    
    /// Main constructor
    /// @param position - position on cylinder 1 and 2, respectively
    /// @param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MLinker(int linkerType, double position1, double position2,
            const vector<double>& coord11, const vector<double>& coord12,
            const vector<double>& coord21, const vector<double>& coord22);
    
    //@{
    /// Getter for constants
    double getStretchingConstant(){return _kStretch;}
    double getEqLength(){return _eqLength;}
    //@}
    
    /// Set parent
    void setLinker(Linker* linker) {_pLinker = linker;}
    /// Get parent 
    Linker* getLinker() {return _pLinker;}
    
    //@{
    /// Length management
    void setLength(double l){_currentLength = l;}
    double getLength() {return _currentLength;}
    //@}
    
private:
    double _eqLength;  ///< Equilibrium length, set at construction
    double _kStretch;  ///< Stretching constant
    
    Linker* _pLinker; ///< Pointer to parent linker
    
    double _currentLength; ///< Current length of the linker
    
};

#endif 
