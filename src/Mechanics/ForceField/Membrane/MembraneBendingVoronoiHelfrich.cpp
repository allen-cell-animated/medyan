#include <array>
#include <vector>

#include "MembraneBendingVoronoiHelfrich.h"

#include "Structure/SurfaceMesh/Vertex.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MembraneBendingVoronoiHelfrich::energy(double area, double curv,
                                              double kBending, double eqCurv){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = curv - eqCurv;
    
    return 2 * kBending * dist * dist * area;
    
}

void MembraneBendingVoronoiHelfrich::forces(
    double* force,
    double area, const mathfunc::Vec3& dArea,
    double curv, const mathfunc::Vec3& dCurv,
    double kBending, double eqCurv
) {
    // F_i = -grad_i U = -4k (H - c0) A (grad_i H) - 2k (H - c0)^2 (grad_i A)
    // A(area), grad_i A(area), H(curv) and grad_i H(curv) are obtained as function parameters

    double dist = curv - eqCurv;
    double coeff1 = -4 * kBending * dist * area;
    double coeff2 = -2 * kBending * dist * dist;

    for(size_t i = 0; i < 3; ++i) force[i] += coeff1 * dCurv[i] + coeff2 * dArea[i];
}
