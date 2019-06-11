#ifndef MEDYAN_Util_Math_ColorMap_Hpp
#define MEDYAN_Util_Math_ColorMap_Hpp

#include <cmath> // abs

#include "Util/Math/Vec.hpp"

// A colormap is a function, which transforms a real number (0.0 - 1.0) to an
// RGB value
namespace colormap {

namespace utility {
    // Utility functions used for colormaps

    // Clamp function (use std::clamp in <algorithm> starting C++17)
    template< typename Float >
    constexpr Float clamp(Float v, Float lo, Float hi) { return (v < lo) ? lo : (v > hi) ? hi : v; }

    template< typename Float >
    constexpr Float clamp01(Float v) { return clamp(v, (Float)0.0, (Float)1.0); }

} // namespace utility

template< typename FloatOut, typename FloatIn >
inline
mathfunc::Vec< 3, FloatOut >
jet(FloatIn v) {
    constexpr auto interpolation = [](FloatIn v) { return clamp01((FloatIn)1.5 - std::abs(2 * v)); };
    return mathfunc::Vec< 3, FloatOut > {
        interpolation(v - 0.5),
        interpolation(v),
        interpolation(v + 0.5)
    };
}

} // namespace colormap

#endif
