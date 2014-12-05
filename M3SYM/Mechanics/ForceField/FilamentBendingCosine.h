//
//  FilamentBendingCosine.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__FilamentBendingCosine__
#define __M3SYM__FilamentBendingCosine__

//FORWARD DECLARATIONS
class Bead;

/// FilamentBendingHarmonic class is a cosine potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentBendingCosine {
    
public:
    double energy(Bead*, Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, Bead*, double, double, double);
    void forces(Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, double, double);
};

#endif /* defined(__M3SYM__FilamentBendingCosine__) */
