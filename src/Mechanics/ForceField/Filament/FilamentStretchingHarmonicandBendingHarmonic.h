
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

#ifndef MEDYAN_FilamentStretchingHarmonicandBendingHarmonic_h
#define MEDYAN_FilamentStretchingHarmonicandBendingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentStretchingHarmonicandBendingHarmonic {

public:
	void energy(floatingpoint *coord, std::size_t nint, int *beadSet,
	            floatingpoint *kstr, floatingpoint *kbend, floatingpoint *eql,
	            floatingpoint *eqt, floatingpoint* totalenergy, const int startID, const
	            int endID, int threadID);
	void energy(floatingpoint *coord, int *beadSetsansbending, floatingpoint *kstrsansbending,
	            floatingpoint *eqlsansbending, floatingpoint* totalenergy, const int startID,
	            const int endID, int threadID);

	void forces(floatingpoint *coord, floatingpoint *f, std::size_t nint, int *beadSet,
	            floatingpoint *kbend, floatingpoint *eqt);
};

#endif
