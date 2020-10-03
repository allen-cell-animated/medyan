#ifndef MEDYAN_Visual_MeshData_hpp
#define MEDYAN_Visual_MeshData_hpp

// This file defines the data structure for mesh data stored in the memory,
// which can then be mapped into GPU memory for drawing.

#include "Visual/FrameData.hpp"

namespace medyan::visual {

//-------------------------------------
// Settings (views)
//-------------------------------------

struct CylindricalObjectDisplaySettings {
    int sides = 10;
};

// The selector is by default choosing all membranes
struct MembraneDisplaySettings {

    bool enabled = true;

    mathfunc::Vec3f colorFixed { 0.4f, 0.6f, 0.95f };
};

//-------------------------------------
// Data structures
//-------------------------------------

struct MeshDataDescriptor {
    int strideSize = 9;

    int positionStart = 0;
    int positionSize = 3;
    int normalStart = 3;
    int normalSize = 3;
    int colorStart = 6;
    int colorSize = 3;
};

struct MeshData {
    MeshDataDescriptor descriptor;
    std::vector< float > data;
};

// preset descriptors
inline constexpr MeshDataDescriptor meshDataDescriptorSurface { 9, 0, 3, 3, 3, 6, 3 };
inline constexpr MeshDataDescriptor meshDataDescriptorLine    { 6, 0, 3, 3, 0, 3, 3 };


//-------------------------------------
// Functions
//-------------------------------------

inline MeshData createMembraneMesh(
    const MembraneFrame& membrane,
    const MembraneDisplaySettings& membraneSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorSurface;

    res.data.reserve(3 * membrane.triangles.size() * res.descriptor.strideSize);
    for(const auto& t : membrane.triangles) {
        const auto& c0 = membrane.vertexCoords[t[0]];
        const auto& c1 = membrane.vertexCoords[t[1]];
        const auto& c2 = membrane.vertexCoords[t[2]];
        const auto un = normalizedVector(cross(c1 - c0, c2 - c0));

        for(size_t i = 0; i < 3; ++i) {
            res.data.push_back(membrane.vertexCoords[t[i]][0]);
            res.data.push_back(membrane.vertexCoords[t[i]][1]);
            res.data.push_back(membrane.vertexCoords[t[i]][2]);
            res.data.push_back(un[0]);
            res.data.push_back(un[1]);
            res.data.push_back(un[2]);
            res.data.push_back(membraneSettings.colorFixed[0]);
            res.data.push_back(membraneSettings.colorFixed[1]);
            res.data.push_back(membraneSettings.colorFixed[2]);
        }
    }

    return res;
}


} // namespace medyan::visual

#endif
