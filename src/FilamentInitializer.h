
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_FilamentInitializer_h
#define MEDYAN_FilamentInitializer_h

#include "common.h"
typedef vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> filamentData;
typedef  tuple< vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> , vector<tuple<string, short, vector<vector<floatingpoint>>>> , vector<tuple<string, short, vector<floatingpoint>>> , vector<vector<floatingpoint>> >  FilamentData;

///FORWARD DECLARATIONS
class Boundary;

/// An interface to initialize an initial configuration of [Filaments](@ref Filament)
/// in the SubSystem.
/*!
 *  FilamentInitiazer class should be inherited to provide an intial scheme for
 *  filling a SubSystem with [Filaments](@ref Filament). The filaments could be 
 *  completely random, directed, a certain length, etc.
 */
class FilamentInitializer {
    
public:
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~FilamentInitializer() noexcept {}
    
    /// Returns a vector of tuples representing the Filament type and beginning and end
    /// coordinates, similar to the structure of manual parsing.
    virtual FilamentData createFilaments(Boundary* b, int numFilaments,
                                                      int filamentType,
                                                      int lenFilaments) = 0;
};

/// An implementation of FilamentInitialzer that creates a completely random
/// Filament distribution within the specified boundary
class RandomFilamentDist : public FilamentInitializer {
    
public:
    FilamentData createFilaments(Boundary* b, int numFilaments,
                                              int filamentType,
                                              int lenFilaments);
};

/// An implementation of FilamentInitialzer that creates a sufficiently spaced
/// network of filaments for investigation of small numbers of filaments.
class ConnectedFilamentDist : public FilamentInitializer {
    
public:
    FilamentData createFilaments(Boundary* b, int numFilaments,
                                 int filamentType,
                                 int lenFilaments);
};

/// An implementation of FilamentInitialzer that creates a random MTOC configuration
class MTOCFilamentDist : public FilamentInitializer {
    
private:
    vector<floatingpoint> _coordMTOC; ///< Coordinates of the MTOC to make filaments around
    floatingpoint _radius;            ///< Radius of MTOC
    
public:
    ///Constructor sets parameters of MTOC
    MTOCFilamentDist(vector<floatingpoint> coord, floatingpoint radius)
        : _coordMTOC(coord), _radius(radius) {}
    
    FilamentData createFilaments(Boundary* b, int numFilaments,
                                              int filamentType,
                                              int lenFilaments);
};


#endif
