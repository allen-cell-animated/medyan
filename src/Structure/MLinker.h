
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_MLinker_h
#define MEDYAN_MLinker_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents the mechanical component of a Linker.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a Linker,
 *  and its associated equilibrium constants. Initial length of a Linker is determined
 *  by the condition of zero initial stress, i.e., it calculated within the constructor
 *  at initiation. A Linker heads positions on a segment (between two consecutive beads 
 *  on a Filament) determined by two numbers (floatingpoint from 0 to 1) position1 and 
 *  position2 held in the parent Linker.
 */
class MLinker {
    
public:
    floatingpoint stretchForce = 0.0; ///< Stretching force of linker at current state
    
    /// Main constructor
    /// @param position - position on cylinder 1 and 2, respectively
    /// @param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MLinker(int linkerType, floatingpoint position1, floatingpoint position2,
            const vector<floatingpoint>& coord11, const vector<floatingpoint>& coord12,
            const vector<floatingpoint>& coord21, const vector<floatingpoint>& coord22);
    
    //@{
    /// Getter for constants
    floatingpoint getStretchingConstant(){return _kStretch;}
    floatingpoint getEqLength(){return _eqLength;}
    //@}
    
    /// Set parent
    void setLinker(Linker* linker) {_pLinker = linker;}
    /// Get parent 
    Linker* getLinker() {return _pLinker;}
    
    //@{
    /// Length management
    void setLength(floatingpoint l){_currentLength = l;}
    floatingpoint getLength() {return _currentLength;}
    //@}
    
private:
    floatingpoint _eqLength;  ///< Equilibrium length, set at construction
    floatingpoint _kStretch;  ///< Stretching constant
    
    Linker* _pLinker; ///< Pointer to parent linker
    
    floatingpoint _currentLength; ///< Current length of the linker
    
};

#endif 
