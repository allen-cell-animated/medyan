#ifndef MEDYAN_Visual_MeshData_hpp
#define MEDYAN_Visual_MeshData_hpp

// This file defines the data structure for mesh data stored in the memory,
// which can then be mapped into GPU memory for drawing.

#include <iterator>
#include <numeric>

#include "Visual/FrameData.hpp"
#include "Visual/Render/PathExtrude.hpp"
#include "Visual/Render/Sphere.hpp"

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

    // path extrude parameters
    float  pathExtrudeRadius = 7.5f;
    int    pathExtrudeSides = 10;

    // bead parameters
    float  beadRadius = 12.0f;
    int    beadLongitudeSegs = 10;
    int    beadLatitudeSegs = 5;

    // display settings
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

// Make mesh data for a single filament and append it to the mesh data.
inline void appendFilamentMeshData(
    MeshData&                      meshData,
    const FilamentFrame&           filament,
    const FilamentDisplaySettings& filamentSettings
) {
    using PM = FilamentDisplaySettings::PathMode;
    using namespace std;

    const int numBeads = filament.coords.size();
    switch(filamentSettings.pathMode) {
        case PM::line:
            if(numBeads > 1) {
                const int numSegments = numBeads - 1;

                for(int i = 0; i < numSegments; ++i) {
                    meshData.data.push_back(filament.coords[i][0]);
                    meshData.data.push_back(filament.coords[i][1]);
                    meshData.data.push_back(filament.coords[i][2]);
                    meshData.data.push_back(filamentSettings.colorFixed[0]);
                    meshData.data.push_back(filamentSettings.colorFixed[1]);
                    meshData.data.push_back(filamentSettings.colorFixed[2]);

                    meshData.data.push_back(filament.coords[i+1][0]);
                    meshData.data.push_back(filament.coords[i+1][1]);
                    meshData.data.push_back(filament.coords[i+1][2]);
                    meshData.data.push_back(filamentSettings.colorFixed[0]);
                    meshData.data.push_back(filamentSettings.colorFixed[1]);
                    meshData.data.push_back(filamentSettings.colorFixed[2]);
                }
            }
            break;

        case PM::extrude:

            {
                vector< int > trivialIndices(numBeads);
                iota(begin(trivialIndices), end(trivialIndices), 0);
                const auto [genVertices, genVertexNormals, genTriInd]
                    = PathExtrude<float> {
                        filamentSettings.pathExtrudeRadius,
                        filamentSettings.pathExtrudeSides
                    }.generate(filament.coords, trivialIndices);

                const int numTriangles = genTriInd.size();
                for(int t = 0; t < numTriangles; ++t) {
                    const auto& triInds = genTriInd[t];
                    const mathfunc::Vec3f coord[] {
                        genVertices[triInds[0]],
                        genVertices[triInds[1]],
                        genVertices[triInds[2]]
                    };
                    const mathfunc::Vec3f un[] {
                        genVertexNormals[triInds[0]],
                        genVertexNormals[triInds[1]],
                        genVertexNormals[triInds[2]]
                    };

                    for(int i = 0; i < 3; ++i) {
                        meshData.data.push_back(coord[i][0]);
                        meshData.data.push_back(coord[i][1]);
                        meshData.data.push_back(coord[i][2]);
                        meshData.data.push_back(un[i][0]);
                        meshData.data.push_back(un[i][1]);
                        meshData.data.push_back(un[i][2]);
                        meshData.data.push_back(filamentSettings.colorFixed[0]);
                        meshData.data.push_back(filamentSettings.colorFixed[1]);
                        meshData.data.push_back(filamentSettings.colorFixed[2]);
                    }
                }
            }
            break;

        case PM::bead:

            {
                const auto sphereGen = SphereUv<float> {
                    filamentSettings.beadRadius,
                    filamentSettings.beadLongitudeSegs,
                    filamentSettings.beadLatitudeSegs
                };
                const auto sphereCache = sphereGen.makeCache();

                for(const auto& bc : filament.coords) {

                    const auto [genVertices, _] = sphereGen.generate(
                        {
                            static_cast<float>(bc[0]),
                            static_cast<float>(bc[1]),
                            static_cast<float>(bc[2])
                        },
                        sphereCache
                    );

                    const int numTriangles = sphereCache.triInd.size();
                    for(int t = 0; t < numTriangles; ++t) {
                        const typename decltype(genVertices)::value_type coord[] {
                            genVertices[sphereCache.triInd[t][0]],
                            genVertices[sphereCache.triInd[t][1]],
                            genVertices[sphereCache.triInd[t][2]]
                        };
                        const auto un = normalizedVector(cross(coord[1] - coord[0], coord[2] - coord[0]));

                        for(int i = 0; i < 3; ++i) {
                            meshData.data.push_back(coord[i][0]);
                            meshData.data.push_back(coord[i][1]);
                            meshData.data.push_back(coord[i][2]);
                            meshData.data.push_back(un[0]);
                            meshData.data.push_back(un[1]);
                            meshData.data.push_back(un[2]);
                            meshData.data.push_back(filamentSettings.colorFixed[0]);
                            meshData.data.push_back(filamentSettings.colorFixed[1]);
                            meshData.data.push_back(filamentSettings.colorFixed[2]);
                        }
                    }
                } // end loop beads in a filament

            }
            break;
    } // end switch

}


inline auto createFilamentMeshData(
    const std::vector<FilamentFrame>& filaments,
    const FilamentDisplaySettings& filamentSettings
) {
    using PM = FilamentDisplaySettings::PathMode;

    MeshData res;

    // reserve space and set descriptor
    switch(filamentSettings.pathMode) {
        case PM::line:
            res.descriptor = meshDataDescriptorLine;
            {
                int numSegments = 0;
                for(auto& filament : filaments) {
                    if(filament.coords.size() >= 2) {
                        numSegments += filament.coords.size() - 1;
                    }
                }
                res.data.reserve(2 * numSegments * res.descriptor.strideSize);
            }
            break;

        case PM::extrude:
            res.descriptor = meshDataDescriptorSurface;
            {
                int numTriangles = 0;
                for(auto& filament : filaments) {
                    numTriangles += PathExtrude<float>{
                        filamentSettings.pathExtrudeRadius,
                        filamentSettings.pathExtrudeSides
                    }.estimateNumTriangles(filament.coords.size());
                }
                res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
            }
            break;

        case PM::bead:
            res.descriptor = meshDataDescriptorSurface;
            {
                const int numTrianglesPerBead = SphereUv<float> {
                    filamentSettings.beadRadius,
                    filamentSettings.beadLongitudeSegs,
                    filamentSettings.beadLatitudeSegs
                }.estimateNumTriangles();

                int numBeads = 0;
                for(auto& filament : filaments) {
                    numBeads += filament.coords.size();
                }

                res.data.reserve(3 * numTrianglesPerBead * numBeads * res.descriptor.strideSize);
            }
            break;
    }

    for(auto& filament : filaments) {
        appendFilamentMeshData(res, filament, filamentSettings);
    }

    return res;
}


} // namespace medyan::visual

#endif
