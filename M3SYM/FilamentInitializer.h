
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

#ifndef M3SYM_FilamentInitializer_h
#define M3SYM_FilamentInitializer_h

#include "common.h"
#include "utility.h"

#include "SubSystem.h"

///FORWARD DECLARATIONS
class Boundary;

/// An interface to initialize an initial configuration of filaments in the system
/*!
 *  FilamentInitiazer class should be inherited to provide an intial scheme for
 *  filling a subsystem with filaments. The filaments could be completely random,
 *  directed, a certain length, etc.
 */

class FilamentInitializer {
    
public:
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~FilamentInitializer() noexcept {}
    
    virtual vector<vector<vector<double>>>
        createFilaments(Boundary* b, int numFilaments, int lenFilaments) = 0;
};

/// An implementation of FilamentInitialzer that creates a completely random
/// filament distribution.
class RandomFilamentDist : public FilamentInitializer {
    
public:
    virtual vector<vector<vector<double>>>
        createFilaments(Boundary* b, int numFilaments, int lenFilaments);
};

#endif
