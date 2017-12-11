#include <array>
#include <vector>

#include "MembraneBendingVoronoiHelfrich.h"

#include "Vertex.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MembraneBendingVoronoiHelfrich::energy(double area, double curv,
                                              double kBending, double eqCurv){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = curv - eqCurv;
    
    return 2 * kBending * dist * dist / area;
    
}

double MembraneBendingVoronoiHelfrich::energy(double areaStretched, double curvStretched,
                                              double kBending, double eqCurv, double d){
    // In fact, d is a dummy variable here, as areaStretched is already dependent on d.

    double distStretched = curvStretched - eqCurv;
    return 2 * kBending * distStretched * distStretched / areaStretched;
}

void MembraneBendingVoronoiHelfrich::forces(Vertex* vCenter, const std::vector<Vertex*>& v,
                                            double area, const std::array<double, 3>& dArea, const std::vector<std::array<double, 3>>& dNeighborArea,
                                            double curv, const std::array<double, 3>& dCurv, const std::vector<std::array<double, 3>>& dNeighborCurv,
                                            double kBending, double eqCurv) {
    // F_i = -grad_i U = -4k (H - c0) A (grad_i H) - 2k (H - c0)^2 (grad_i A)
    // A(area), grad_i A(area), H(curv) and grad_i H(curv) are obtained as function parameters

    size_t n = v.size();

    double dist = curv - eqCurv;
    double coeff1 = -4 * kBending * dist * area;
    double coeff2 = -2 * kBending * dist * dist;

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        vCenter->force[coordIdx] += coeff1 * dArea[coordIdx] + coeff2 * dCurv[coordIdx];
        for(size_t nIdx = 0; nIdx < n; ++nIdx) {
            v[nIdx]->force[coordIdx] += coeff1 * dNeighborArea[nIdx][coordIdx] + coeff2 * dNeighborCurv[nIdx][coordIdx];
        }
    }
}

void MembraneBendingVoronoiHelfrich::forcesAux(Vertex* vCenter, const std::vector<Vertex*>& v,
                                               double area, const std::array<double, 3>& dArea, const std::vector<std::array<double, 3>>& dNeighborArea,
                                               double curv, const std::array<double, 3>& dCurv, const std::vector<std::array<double, 3>>& dNeighborCurv,
                                               double kBending, double eqCurv) {
    // Same as force calculation

    size_t n = v.size();

    double dist = curv - eqCurv;
    double coeff1 = -4 * kBending * dist * area;
    double coeff2 = -2 * kBending * dist * dist;

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        vCenter->forceAux[coordIdx] += coeff1 * dArea[coordIdx] + coeff2 * dCurv[coordIdx];
        for(size_t nIdx = 0; nIdx < n; ++nIdx) {
            v[nIdx]->forceAux[coordIdx] += coeff1 * dNeighborArea[nIdx][coordIdx] + coeff2 * dNeighborCurv[nIdx][coordIdx];
        }
    }
}
