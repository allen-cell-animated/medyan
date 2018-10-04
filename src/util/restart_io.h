#ifndef MEDYAN_UTIL_RESTART_IO_H
#define MEDYAN_UTIL_RESTART_IO_H

#include <algorithm>
#include <fstream>
#include <type_traits>
#include <vector>

#include "util/restart_data.h"

// Tasks for restart io:
//
// Reading:
//   - Read restart files and store them in restart data and/or SysParams
//
// Writing:
//   - Write restart files from restart data and/or SysParams

namespace restart {

// "unsigned char"s should be replaced by "std::byte"s starting C++17
using vec_byte = std::vector<unsigned char>;

namespace helper {

template<
    typename T,
    typename = typename std::enable_if< std::is_standard_layout< T >::value >::type
> struct StdLayoutIntepreter {
    vec_byte raw(const T& x)const {
        unsigned char const * pt = reinterpret_cast<unsigned char const*>(&x);
        return vec_byte(pt, pt + sizeof(T));
    }
    void store(T& x, const vec_byte& b)const {
        // Must satisfy b.size() == sizeof(T). Otherwise UB.
        unsigned char *pt = reinterpret_cast<unsigned char*>(&x);
        std::copy(b.begin(), b.end(), pt);
    }
};

} // namespace helper

struct Layer {
    std::fstream& fs;

    RestartData& rd;

    vec_byte header;
    // body could be Layer, a vector of Layers or pure bytes

    // Read from restart file (from current ofstream) and write information to system
    virtual void read() = 0;

    // Read from system and write information to restart file (from current ifstream)
    virtual void write() = 0;
};

struct LayerElement : public Layer {
};
struct LayerFilament : public Layer {
    RestartData::FilamentData& thisFilamentData;

    virtual void read() override;
    virtual void write() override;
};

} // namespace restart

#endif // include guard
