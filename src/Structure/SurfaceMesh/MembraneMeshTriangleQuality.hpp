#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshTriangleQuality_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshTriangleQuality_hpp

#include <limits>

#include "MathFunctions.h"

enum class TriangleQualityCriteria {
    RadiusRatio     // Circumradius / (2 * Inradius), range [1, inf)
};
template< TriangleQualityCriteria > struct TriangleQuality;
template<> struct TriangleQuality< TriangleQualityCriteria::RadiusRatio > {
    static constexpr double best = 1.0;
    static constexpr double worst = std::numeric_limits<double>::infinity();
    static constexpr bool better(double q1, double q2) { return q1 < q2; }
    static constexpr auto betterOne(double q1, double q2) { return better(q1, q2) ? q1 : q2; }
    static constexpr bool worse(double q1, double q2) { return q1 > q2; }
    static constexpr auto worseOne(double q1, double q2) { return worse(q1, q2) ? q1 : q2; }
    static constexpr auto improvement(double q0, double q1) { return q0 / q1; }

    template< typename VecType >
    auto operator()(const VecType& v0, const VecType& v1, const VecType& v2) const {
        using namespace mathfunc;
        const auto d0 = distance(v1, v2);
        const auto d1 = distance(v2, v0);
        const auto d2 = distance(v0, v1);
        return operator()(d0, d1, d2);
    }

    auto operator()(double d0, double d1, double d2) const {
        const auto p = 0.5 * (d0 + d1 + d2);
        // Note that the abs is needed to avoid extreme cases where the result is negative.
        return std::abs(d0 * d1 * d2 / (8 * (p - d0) * (p - d1) * (p - d2)));
    }
};

#endif
