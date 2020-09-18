#ifndef MEDYAN_Visual_FrameData_hpp
#define MEDYAN_Visual_FrameData_hpp

#include <array>
#include <cstddef>
#include <vector>

#include "Util/Math/Vec.hpp"

namespace medyan::visual {


//-------------------------------------
// Data structures
//-------------------------------------

struct MembraneFrame {
    // meta data
    int id = 0;

    // triangular mesh data
    std::vector< std::array< int, 3 >> triangles;
    std::vector< mathfunc::Vec3d > vertexCoords;
};


struct DisplayFrame {
    // meta data

    // frame data
    MembraneFrame membrane;
};


//-------------------------------------
// Functions
//-------------------------------------

inline DisplayFrame readFrameDataFromOutput() {
    DisplayFrame res;

    // Not implemented yet

    return res;
}


} // namespace medyan::visual

#endif
