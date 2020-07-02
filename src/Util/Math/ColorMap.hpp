#ifndef MEDYAN_Util_Math_ColorMap_Hpp
#define MEDYAN_Util_Math_ColorMap_Hpp

#include <algorithm>
#include <array>
#include <cmath> // abs

#include "Util/Math/Vec.hpp"

// A colormap is a function, which transforms a real number in a range to an
// RGB value
namespace colormap {

template< typename Float >
using ColorRgb = mathfunc::Vec< 3, Float >;

enum Extend {
    clamp,
    cycle,
    bounce
};

// Color constants
template< typename Float >
inline constexpr ColorRgb< Float > black { (Float)0.0, (Float)0.0, (Float)0.0 };
inline constexpr auto blackD = black< double >;
inline constexpr auto blackF = black< float >;
template< typename Float >
inline constexpr ColorRgb< Float > white { (Float)1.0, (Float)1.0, (Float)1.0 };
inline constexpr auto whiteD = white< double >;
inline constexpr auto whiteF = white< float >;

template< typename Float >
inline constexpr ColorRgb< Float > red { (Float)1.0, (Float)0.0, (Float)0.0 };
inline constexpr auto redD = red< double >;
inline constexpr auto redF = red< float >;
template< typename Float >
inline constexpr ColorRgb< Float > green { (Float)0.0, (Float)1.0, (Float)0.0 };
inline constexpr auto greenD = green< double >;
inline constexpr auto greenF = green< float >;
template< typename Float >
inline constexpr ColorRgb< Float > blue { (Float)0.0, (Float)0.0, (Float)1.0 };
inline constexpr auto blueD = blue< double >;
inline constexpr auto blueF = blue< float >;


// Color range extension
template< typename Float >
constexpr Float clamp01(Float v) { return std::clamp(v, (Float)0.0, (Float)1.0); }

// Not constexpr because cmath functions are not constexpr
template< typename Float >
inline Float cycle01(Float v) {
    return v - std::floor(v);
}

template< typename Float >
inline Float bounce01(Float v) {
    const auto v1 = v - 2 * std::floor(v / 2);
    return v > 1 ? 2 - v : v;
}

template< typename Float >
inline Float extend01(Float v, Extend ex) {
    switch(ex) {
    case Extend::clamp:
        return clamp01(v);
    case Extend::cycle:
        return cycle01(v);
    case Extend::bounce:
        return bounce01(v);
    }
}

// Color interpolation
// Precondition: v must be in range [0, 1].
template< typename FloatVal, size_t colorDim, typename FloatColor, size_t listSize >
constexpr
mathfunc::Vec< colorDim, FloatColor >
interpolate(FloatVal v, std::array< mathfunc::Vec< colorDim, FloatColor >, listSize > interpList) {
    static_assert(listSize > 1, "Must have at least 2 elements for interpolation");
    const FloatVal interval = static_cast< FloatVal >(1.0) / (listSize - 1);
    const FloatVal normalizedPos = v / interval;
    const size_t intervalIndex = static_cast< size_t >(normalizedPos);
    if(intervalIndex >= listSize - 1) return interpList[listSize - 1];

    return interpList[intervalIndex    ] * static_cast< FloatColor >(intervalIndex + 1 - normalizedPos)
        +  interpList[intervalIndex + 1] * static_cast< FloatColor >(normalizedPos - intervalIndex);
}


// Define color maps
//-----------------------------------------------------------------------------

template< typename FloatColor >
struct Jet {
    static constexpr std::array< ColorRgb< FloatColor >, 9 > interpList = {{
        { 0.0, 0.0, 0.5 },
        { 0.0, 0.0, 1.0 },
        { 0.0, 0.5, 1.0 },
        { 0.0, 1.0, 1.0 },
        { 0.5, 1.0, 0.5 },
        { 1.0, 1.0, 0.0 },
        { 1.0, 0.5, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.5, 0.0, 0.0 }
    }};
};
using JetD = Jet< double >;
using JetF = Jet< float >;

template< typename FloatColor >
struct Bwr {
    static constexpr std::array< ColorRgb< FloatColor >, 3 > interpList = {{
        { 0.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 1.0, 0.0, 0.0 }
    }};
};
using BwrD = Bwr< double >;
using BwrF = Bwr< float >;

// Get the color
template< typename Map, typename FloatVal >
constexpr auto color(Map, FloatVal v, Extend ex = Extend::clamp) {
    return interpolate(extend01(v), Map::interpList);
}

} // namespace colormap

#endif
