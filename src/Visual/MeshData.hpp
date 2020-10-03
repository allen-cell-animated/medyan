#ifndef MEDYAN_Visual_MeshData_hpp
#define MEDYAN_Visual_MeshData_hpp

// This file defines the data structure for mesh data stored in the memory,
// which can then be mapped into GPU memory for drawing.

#include "Visual/FrameData.hpp"

namespace medyan::visual {

//-------------------------------------
// Settings (for display)
//-------------------------------------

struct SurfaceDisplaySettings {
    enum class PolygonMode { wireframe, fill };

    bool enabled = true;
    PolygonMode polygonMode = PolygonMode::wireframe;
};

struct LineDisplaySettings {
    bool enabled = true;
};


//-------------------------------------
// Settings (for making mesh data)
//-------------------------------------

struct MembraneDisplaySettings {
    mathfunc::Vec3f colorFixed { 0.4f, 0.6f, 0.95f };

    SurfaceDisplaySettings surface;
};

struct FilamentDisplaySettings {
    enum class PathMode { line, extrude, bead };

    mathfunc::Vec3f colorFixed { 0.95f, 0.1f, 0.15f };

    PathMode pathMode = PathMode::extrude;
    SurfaceDisplaySettings surface;  // used in "extrude" or "bead" path mode
    LineDisplaySettings    line;     // used in "line" path mode

    bool isEnabled() const {
        return pathMode == PathMode::line ?
            line.enabled :
            surface.enabled;
    }
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

inline void appendMembraneMeshData(
    MeshData&                      meshData,
    const MembraneFrame&           membrane,
    const MembraneDisplaySettings& membraneSettings
) {

    for(const auto& t : membrane.triangles) {
        const auto& c0 = membrane.vertexCoords[t[0]];
        const auto& c1 = membrane.vertexCoords[t[1]];
        const auto& c2 = membrane.vertexCoords[t[2]];
        const auto un = normalizedVector(cross(c1 - c0, c2 - c0));

        for(size_t i = 0; i < 3; ++i) {
            meshData.data.push_back(membrane.vertexCoords[t[i]][0]);
            meshData.data.push_back(membrane.vertexCoords[t[i]][1]);
            meshData.data.push_back(membrane.vertexCoords[t[i]][2]);
            meshData.data.push_back(un[0]);
            meshData.data.push_back(un[1]);
            meshData.data.push_back(un[2]);
            meshData.data.push_back(membraneSettings.colorFixed[0]);
            meshData.data.push_back(membraneSettings.colorFixed[1]);
            meshData.data.push_back(membraneSettings.colorFixed[2]);
        }
    }

}

inline auto createMembraneMeshData(
    const MembraneFrame& membrane,
    const MembraneDisplaySettings& membraneSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorSurface;

    res.data.reserve(3 * membrane.triangles.size() * res.descriptor.strideSize);
    appendMembraneMeshData(res, membrane, membraneSettings);

    return res;
}

inline auto createMembraneMeshData(
    const std::vector<MembraneFrame>& membranes,
    const MembraneDisplaySettings& membraneSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorSurface;

    int numTriangles = 0;
    for(auto& membrane : membranes) numTriangles += membrane.triangles.size();
    res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
    for(auto& membrane : membranes) {
        appendMembraneMeshData(res, membrane, membraneSettings);
    }

    return res;
}

inline void appendFilamentMeshData(
    MeshData&                      meshData,
    const FilamentFrame&           filament,
    const FilamentDisplaySettings& filamentSettings
) {

}


} // namespace medyan::visual

#endif
