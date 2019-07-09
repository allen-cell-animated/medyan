#ifndef MEDYAN_TESTS_Mechanics_ForceField_TestFFCommon_Hpp
#define MEDYAN_TESTS_Mechanics_ForceField_TestFFCommon_Hpp

#include <cmath>

#include "Rand.h"
#include "Util/Math/Vec.hpp"

// Provides common methods for force field testing purposes
namespace test_ff_common {

// Fill vectors with random numbers
template< typename ForwardIt, typename Float >
inline void fillNormalRand(ForwardIt first, ForwardIt last, Float mean, Float stddev) {
    std::normal_distribution< Float > nd(mean, stddev);
    std::generate(first, last, [&nd] { return nd(Rand::eng); });
}

template< size_t dim, typename Float, typename Container >
inline void fillNormalRand(mathfunc::VecArray< dim, Float, Container >& v, Float mean, Float stddev) {
    fillNormalRand(v.value.begin(), v.value.end(), mean, stddev);
}

// Vector increment
template< size_t dim, typename Float, typename Container >
inline void vecInc(mathfunc::VecArray< dim, Float, Container >& v1, const mathfunc::VecArray< dim, Float, Container >& v2) {
    v1 += v2;
}

// Vector decrement
template< size_t dim, typename Float, typename Container >
inline void vecDec(mathfunc::VecArray< dim, Float, Container >& v1, const mathfunc::VecArray< dim, Float, Container >& v2) {
    v1 -= v2;
}

// Dot product
template< size_t dim, typename Float, typename Container >
inline Float vecDot(mathfunc::VecArray< dim, Float, Container >& v1, const mathfunc::VecArray< dim, Float, Container >& v2) {
    return mathfunc::dot(v1, v2);
}

// Compare equal
template< typename Float >
inline bool equalRelEps(Float a, Float b, Float relEps) {
    return std::fabs(a - b) <= std::min(std::fabs(a), std::fabs(b)) * relEps;
}

// Energy - force consistency test
template< typename CoordContainer, typename Float >
struct EnergyForceConsistencyReport {
    bool passed;
    CoordContainer dc;
    Float deActual;
    Float deExpected;
};

template<
    typename CoordContainer,
    typename Float,
    typename FuncEnergy,
    typename FuncForce
> inline auto testEnergyForceConsistency(
    const CoordContainer& c0,
    FuncEnergy&& calcEnergy, FuncForce&& calcForce,
    Float moveMag, Float energyRelEps
) {
    EnergyForceConsistencyReport< CoordContainer, Float > res;

    res.dc.resize(c0.size());
    fillNormalRand(res.dc, (Float)0.0, moveMag);
    CoordContainer c1 = c0; vecDec(c1, res.dc);
    CoordContainer c2 = c0; vecInc(c2, res.dc);

    CoordContainer f;
    f.resize(c0.size());

    // Actual change in energy
    const Float e1 = calcEnergy(c1);
    const Float e2 = calcEnergy(c2);
    res.deActual = e2 - e1;

    // Expected change in energy
    calcForce(c0, f);
    res.deExpected = -2 * vecDot(res.dc, f);

    // Judge equality
    res.passed = equalRelEps(res.deActual, res.deExpected, energyRelEps);

    return res;
}

} // namespace test_ff_common

#endif
