#ifndef MEDYAN_UTIL_MATH_VEC_HPP
#define MEDYAN_UTIL_MATH_VEC_HPP

#include <array>
#include <ostream>

namespace mathfunc {

// A simple coordinate type that makes operations easier
// Can be replaced by linalg/tensor libraries in the future
template< size_t dim, typename Float = double > struct Vec {

    using storage_type = std::array< Float, dim >;
    using size_type = typename storage_type::size_type;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    storage_type value;

    constexpr iterator       begin()       noexcept { return value.begin(); }
    constexpr const_iterator begin() const noexcept { return value.begin(); }
    constexpr iterator       end()       noexcept { return value.end(); }
    constexpr const_iterator end() const noexcept { return value.end(); }

    constexpr       Float& operator[](size_type pos)       { return value[pos]; }
    constexpr const Float& operator[](size_type pos) const { return value[pos]; }
};

// Frequently used type alias
using Vec3 = Vec<3>;

// Formatting
// undefined behavior if dim is 0
template< size_t dim, typename Float > inline
std::ostream& operator<<(std::ostream& os, const Vec<dim, Float>& v) {
    os << '(';
    for(size_t i = 0; i < dim - 1; ++i) os << v[i] << ", ";
    os << v.v[dim - 1] << ')';
    return os;
}

} // namespace mathfunc

#endif
