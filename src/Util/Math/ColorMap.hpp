#ifndef MEDYAN_Util_Math_ColorMap_Hpp
#define MEDYAN_Util_Math_ColorMap_Hpp

#include <algorithm>
#include <array>
#include <cmath> // abs

#include "Util/Math/Vec.hpp"

// A colormap is a function, which transforms a real number in a range to an
// RGB value
namespace colormap {

namespace utility {
    // Utility functions used for colormaps

    template< typename Float >
    constexpr Float clamp01(Float v) { return std::clamp(v, (Float)0.0, (Float)1.0); }

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

} // namespace utility

template< typename FloatColor >
struct Jet {
    static constexpr std::array< mathfunc::Vec< 3, FloatColor >, 9 > interpList = {{
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

    template< typename FloatIn >
    constexpr auto operator()(FloatIn v) const {
        return utility::interpolate(v, interpList);
    }
};

} // namespace colormap

#endif
