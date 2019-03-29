#ifndef MEDYAN_MembraneBendingVoronoiHelfrich_h
#define MEDYAN_MembraneBendingVoronoiHelfrich_h

#include <array>
#include <vector>

#include "common.h"
#include "MathFunctions.h"

//FORWARD DECLARATIONS
class Vertex;

/// A harmonic potential used by the MembraneStretching
class MembraneBendingVoronoiHelfrich {
    
public:
    double energy(double, double, double, double);
    
    void forces(double* force,
        double area, const mathfunc::Vec3& dArea,
        double curv, const mathfunc::Vec3& dCurv,
        double, double);
    void forcesAux(Vertex*,
        double area, const mathfunc::Vec3& dArea,
        double curv, const mathfunc::Vec3& dCurv,
        double, double);
};

#endif
