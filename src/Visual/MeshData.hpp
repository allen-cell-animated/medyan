#ifndef MEDYAN_Visual_MeshData_hpp
#define MEDYAN_Visual_MeshData_hpp

// This file defines the data structure for mesh data stored in the memory,
// which can then be mapped into GPU memory for drawing.

#include <iterator>
#include <numeric>
#include <tuple>

#include "Visual/Common.hpp"
#include "Visual/FrameData.hpp"
#include "Visual/Render/PathExtrude.hpp"
#include "Visual/Render/Sphere.hpp"
#include "Visual/Shader.hpp"

namespace medyan::visual {

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
// preset descriptors
inline constexpr MeshDataDescriptor meshDataDescriptorSurface { 9, 0, 3, 3, 3, 6, 3 };
inline constexpr MeshDataDescriptor meshDataDescriptorLine    { 6, 0, 3, 3, 0, 3, 3 };

class GlVertexBufferManager {
private:
    // vao, vbo, ebo
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    // GLuint ebo_;

    bool bufferMade_ = false;

public:

    auto vao() const { return vao_; }
    auto vbo() const { return vbo_; }

    // ctor and dtor
    GlVertexBufferManager() = default;
    GlVertexBufferManager(MeshDataDescriptor desc) {
        makeBuffer(desc);
    }
    GlVertexBufferManager(GlVertexBufferManager&& rhs) {
        swap(*this, rhs);
    }

    ~GlVertexBufferManager() {
        deleteBuffer();
    }

    GlVertexBufferManager& operator=(GlVertexBufferManager&& rhs) {
        if(this != &rhs) {
            // clear this
            deleteBuffer();

            swap(*this, rhs);
        }
    }

    friend void swap(GlVertexBufferManager& lhs, GlVertexBufferManager& rhs) {
        std::swap(lhs.vao_,        rhs.vao_);
        std::swap(lhs.vbo_,        rhs.vbo_);
        std::swap(lhs.bufferMade_, rhs.bufferMade_);
    }

    void makeBuffer(MeshDataDescriptor desc) {
        if(bufferMade_) deleteBufferImpl_();

        makeBufferImpl_(desc);
        bufferMade_ = true;
    }
    void deleteBuffer() {
        if(bufferMade_) {
            deleteBufferImpl_();
        }
        bufferMade_ = false;
    }

private:
    // Note: this function does not delete previously built buffer
    void makeBufferImpl_(MeshDataDescriptor desc) {
        glGenBuffers(1, &vbo_);
        // glGenBuffers(1, &ebo);
        glGenVertexArrays(1, &vao_);

        glBindVertexArray(vao_);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);

        // Vertex attribute
        //---------------------------------------------------------------------
        GLuint idx = 0;
        const auto addBlock = [&](int blockOffset, int blockSize, int strideSize) {
            if(blockSize > 0) {
                glVertexAttribPointer(
                    idx,
                    blockSize,
                    GL_FLOAT,
                    GL_FALSE,
                    strideSize * sizeof(float),
                    static_cast<const char*>(0) + sizeof(float) * blockOffset
                );
                glEnableVertexAttribArray(idx);
                ++idx;
            }
        };
        // Position
        addBlock(desc.positionStart, desc.positionSize, desc.strideSize);
        addBlock(desc.normalStart,   desc.normalSize,   desc.strideSize);
        addBlock(desc.colorStart,    desc.colorSize,    desc.strideSize);

        // temporarily retarget
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    // Note: this function does not check for creation
    void deleteBufferImpl_() {
        glDeleteVertexArrays(1, &vao_);
        glDeleteBuffers(1, &vbo_);
        // glDeleteBuffers(1, &ebo_);
    }
};


struct MeshData {
    MeshDataDescriptor descriptor;
    std::vector< float > data;

    // Auxiliary state to indicate that the data is updated.
    // Can be set to false by the downstream pipeline.
    bool updated = true;
};




//-------------------------------------
// Settings (for display)
//-------------------------------------

enum class DisplayGeometryType { surface, line };

struct SurfaceDisplaySettings {
    enum class PolygonMode { wireframe, fill };

    bool            enabled = true;
    PolygonMode     polygonMode = PolygonMode::fill;

    mathfunc::Vec3f colorSpecular { 0.5, 0.5, 0.5 };
    float           colorShininess = 32.0f;

    GlVertexBufferManager vertexBufferManager;

    SurfaceDisplaySettings() : vertexBufferManager(meshDataDescriptorSurface) {}
};

struct LineDisplaySettings {
    bool enabled = true;

    GlVertexBufferManager vertexBufferManager;

    LineDisplaySettings() : vertexBufferManager(meshDataDescriptorLine) {}
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

};

struct LinkerDisplaySettings {
    enum class PathMode { line, extrude };

    mathfunc::Vec3f colorFixed { 0.1f, 0.9f, 0.0f };

    PathMode pathMode = PathMode::extrude;

    // path extrude parameters
    float  pathExtrudeRadius = 7.5f;
    int    pathExtrudeSides = 10;

    // display settings
    SurfaceDisplaySettings surface;  // used in "extrude" or "bead" path mode
    LineDisplaySettings    line;     // used in "line" path mode

};

struct AuxLineDisplaySettings {
    using Flag = std::uint_fast8_t;
    inline static constexpr Flag targetCompartmentBorder = 1 << 0;
    inline static constexpr Flag targetCompartmentAll    = 1 << 1;

    Flag flag = 0;

    mathfunc::Vec3f colorFixed { 1.0f, 1.0f, 1.0f };

    LineDisplaySettings line;
};



//-------------------------------------
// Functions (mesh data creation)
//-------------------------------------

// Note that the following functions may be categorized as the following:
//   - The actual function appending data of each element to the mesh data.
//   - The function appending data of all selected elements of a frame to existing mesh data.
//   - The function creating a new mesh data solely for the data of all selected elements of a frame.

// Returns how many numbers added
inline auto appendMembraneMeshData(
    MeshData&                      meshData,
    const MembraneFrame&           membrane,
    const MembraneDisplaySettings& membraneSettings
) {

    const int sizePrev = meshData.data.size();

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

    const int sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}

// Note:
//   - This function does not set descriptor.
inline auto appendMembraneMeshData(
    MeshData&                         meshData,
    std::vector<const MembraneFrame*> pMembranes,
    const MembraneDisplaySettings&    membraneSettings
) {
    int numAdded = 0;
    for(auto pm : pMembranes) {
        numAdded += appendMembraneMeshData(meshData, *pm, membraneSettings);
    }

    return numAdded;
}

// Generate membrane mesh with a selected range of membranes.
//
// The selected membranes should be replaced with a range object with C++20.
inline auto createMembraneMeshData(
    std::vector<const MembraneFrame*> pMembranes,
    const MembraneDisplaySettings&    membraneSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorSurface;

    int numTriangles = 0;
    for(auto pm : pMembranes) numTriangles += pm->triangles.size();
    res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
    for(auto pm : pMembranes) {
        appendMembraneMeshData(res, *pm, membraneSettings);
    }

    return res;
}

// Make mesh data for a single filament and append it to the mesh data.
//
// Returns how many numbers added.
inline auto appendFilamentMeshData(
    MeshData&                      meshData,
    const FilamentFrame&           filament,
    const FilamentDisplaySettings& filamentSettings
) {
    using PM = FilamentDisplaySettings::PathMode;
    using namespace std;

    const int sizePrev = meshData.data.size();
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

    const int sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}

// Note:
//   - This function does not set descriptor.
inline auto appendFilamentMeshData(
    MeshData&                         meshData,
    std::vector<const FilamentFrame*> pFilaments,
    const FilamentDisplaySettings&    filamentSettings
) {
    int numAdded = 0;
    for(auto pf : pFilaments) {
        numAdded += appendFilamentMeshData(meshData, *pf, filamentSettings);
    }

    return numAdded;
}

// Generate filament mesh with selected filaments
//
// The selected membranes should be replaced with a range object with C++20.
inline auto createFilamentMeshData(
    std::vector<const FilamentFrame*> pFilaments,
    const FilamentDisplaySettings&    filamentSettings
) {
    using PM = FilamentDisplaySettings::PathMode;

    MeshData res;

    // reserve space and set descriptor
    switch(filamentSettings.pathMode) {
        case PM::line:
            res.descriptor = meshDataDescriptorLine;
            {
                int numSegments = 0;
                for(auto pf : pFilaments) {
                    if(pf->coords.size() >= 2) {
                        numSegments += pf->coords.size() - 1;
                    }
                }
                res.data.reserve(2 * numSegments * res.descriptor.strideSize);
            }
            break;

        case PM::extrude:
            res.descriptor = meshDataDescriptorSurface;
            {
                int numTriangles = 0;
                for(auto pf : pFilaments) {
                    numTriangles += PathExtrude<float>{
                        filamentSettings.pathExtrudeRadius,
                        filamentSettings.pathExtrudeSides
                    }.estimateNumTriangles(pf->coords.size());
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
                for(auto pf : pFilaments) {
                    numBeads += pf->coords.size();
                }

                res.data.reserve(3 * numTrianglesPerBead * numBeads * res.descriptor.strideSize);
            }
            break;
    }

    for(auto pf : pFilaments) {
        appendFilamentMeshData(res, *pf, filamentSettings);
    }

    return res;
}


// Make mesh data for a single filament and append it to the mesh data.
//
// Returns how many numbers added.
inline auto appendLinkerMeshData(
    MeshData&                    meshData,
    const LinkerFrame&           linker,
    const LinkerDisplaySettings& linkerSettings
) {
    using PM = LinkerDisplaySettings::PathMode;
    using namespace std;

    const int sizePrev = meshData.data.size();
    switch(linkerSettings.pathMode) {
        case PM::line:
            {
                meshData.data.push_back(linker.coords[0][0]);
                meshData.data.push_back(linker.coords[0][1]);
                meshData.data.push_back(linker.coords[0][2]);
                meshData.data.push_back(linkerSettings.colorFixed[0]);
                meshData.data.push_back(linkerSettings.colorFixed[1]);
                meshData.data.push_back(linkerSettings.colorFixed[2]);

                meshData.data.push_back(linker.coords[1][0]);
                meshData.data.push_back(linker.coords[1][1]);
                meshData.data.push_back(linker.coords[1][2]);
                meshData.data.push_back(linkerSettings.colorFixed[0]);
                meshData.data.push_back(linkerSettings.colorFixed[1]);
                meshData.data.push_back(linkerSettings.colorFixed[2]);
            }
            break;

        case PM::extrude:

            {
                vector< int > trivialIndices { 0, 1 };
                const auto [genVertices, genVertexNormals, genTriInd]
                    = PathExtrude<float> {
                        linkerSettings.pathExtrudeRadius,
                        linkerSettings.pathExtrudeSides
                    }.generate(linker.coords, trivialIndices);

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
                        meshData.data.push_back(linkerSettings.colorFixed[0]);
                        meshData.data.push_back(linkerSettings.colorFixed[1]);
                        meshData.data.push_back(linkerSettings.colorFixed[2]);
                    }
                }
            }
            break;

    } // end switch

    const int sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}


// Note:
//   - This function does not set descriptor.
inline auto appendLinkerMeshData(
    MeshData&                       meshData,
    std::vector<const LinkerFrame*> pLinkers,
    const LinkerDisplaySettings&    linkerSettings
) {
    int numAdded = 0;
    for(auto pl : pLinkers) {
        numAdded += appendLinkerMeshData(meshData, *pl, linkerSettings);
    }

    return numAdded;
}


// Create linkers given the range of linkers
//
// The selected membranes should be replaced with a range object with C++20.
inline auto createLinkerMeshData(
    std::vector<const LinkerFrame*> pLinkers,
    const LinkerDisplaySettings&    linkerSettings
) {
    using PM = LinkerDisplaySettings::PathMode;

    MeshData res;

    // reserve space and set descriptor
    switch(linkerSettings.pathMode) {
        case PM::line:
            res.descriptor = meshDataDescriptorLine;
            {
                const int numSegments = pLinkers.size();
                res.data.reserve(2 * numSegments * res.descriptor.strideSize);
            }
            break;

        case PM::extrude:
            res.descriptor = meshDataDescriptorSurface;
            {
                const int numTriangles = pLinkers.size() * PathExtrude<float>{
                    linkerSettings.pathExtrudeRadius,
                    linkerSettings.pathExtrudeSides
                }.estimateNumTriangles(2);
                res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
            }
            break;
    }

    for(auto pl : pLinkers) {
        appendLinkerMeshData(res, *pl, linkerSettings);
    }

    return res;
}


// Returns how many numbers added.
inline auto appendAuxLineMeshData(
    MeshData&                       meshData,
    const DisplayFrame&             frameData,
    const AuxLineDisplaySettings&   auxLineSettings
) {
    const int sizePrev = meshData.data.size();

    if(frameData.compartmentInfo.has_value()) {

        // Function to generate auxiliary lines
        //
        // Inputs:
        //   - fixedAxis: the axis parallel to the line to be drawn. range 0, 1, 2.
        //   - dax1:      the number of compartments interval along the next axis after fixedAxis.
        //   - dax2:      the number of compartments interval along the 2nd next axis after fixedAxis.
        //
        // Notes:
        //   - The next axis index are found using the rotation 0 -> 1 -> 2 -> 0.
        const auto genLineArray = [&](size_t fixedAxis, size_t dax1, size_t dax2) {
            const int ax1 = (fixedAxis + 1) % 3;
            const int ax2 = (fixedAxis + 2) % 3;
            for(size_t x1 = 0; x1 <= frameData.compartmentInfo->number[ax1]; x1 += dax1) {
                for(size_t x2 = 0; x2 <= frameData.compartmentInfo->number[ax2]; x2 += dax2) {
                    const auto v1 = frameData.compartmentInfo->size[ax1] * x1;
                    const auto v2 = frameData.compartmentInfo->size[ax2] * x2;
                    Vec3f coord0, coord1;
                    coord0[fixedAxis] = 0;
                    coord0[ax1] = v1;
                    coord0[ax2] = v2;
                    coord0 += frameData.compartmentInfo->offset;
                    coord1[fixedAxis] = frameData.compartmentInfo->size[fixedAxis] * frameData.compartmentInfo->number[fixedAxis];
                    coord1[ax1] = v1;
                    coord1[ax2] = v2;
                    coord1 += frameData.compartmentInfo->offset;

                    meshData.data.push_back(coord0[0]);
                    meshData.data.push_back(coord0[1]);
                    meshData.data.push_back(coord0[2]);
                    meshData.data.push_back(auxLineSettings.colorFixed[0]);
                    meshData.data.push_back(auxLineSettings.colorFixed[1]);
                    meshData.data.push_back(auxLineSettings.colorFixed[2]);
                    meshData.data.push_back(coord1[0]);
                    meshData.data.push_back(coord1[1]);
                    meshData.data.push_back(coord1[2]);
                    meshData.data.push_back(auxLineSettings.colorFixed[0]);
                    meshData.data.push_back(auxLineSettings.colorFixed[1]);
                    meshData.data.push_back(auxLineSettings.colorFixed[2]);
                }
            }
        };

        if(auxLineSettings.flag & AuxLineDisplaySettings::targetCompartmentBorder) {
            genLineArray(0, frameData.compartmentInfo->number[1], frameData.compartmentInfo->number[2]);
            genLineArray(1, frameData.compartmentInfo->number[2], frameData.compartmentInfo->number[0]);
            genLineArray(2, frameData.compartmentInfo->number[0], frameData.compartmentInfo->number[1]);
        }
        if(auxLineSettings.flag & AuxLineDisplaySettings::targetCompartmentAll) {
            genLineArray(0, 1, 1);
            genLineArray(1, 1, 1);
            genLineArray(2, 1, 1);
        }
    }

    const int sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}
// Auxiliary lines do not need a selector.
// The actual "selection" happens inside a mesh.

inline auto createAuxLineMeshData(
    const DisplayFrame&             frameData,
    const AuxLineDisplaySettings&   auxLineSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorLine;
    appendAuxLineMeshData(res, frameData, auxLineSettings);

    return res;
}

//-------------------------------------
// Functions (mesh data display)
//-------------------------------------


inline auto convertToGlm(const mathfunc::Vec3f& vec) {
    return glm::vec3( vec[0], vec[1], vec[2] );
}


// Auxiliary functions to get polygon mode
inline auto polygonModeGl(const SurfaceDisplaySettings& settings) {
    return settings.polygonMode == SurfaceDisplaySettings::PolygonMode::wireframe ? GL_LINE : GL_FILL;
}
inline auto polygonModeGl(const LineDisplaySettings& settings) {
    return GL_LINE;
}

// Auxiliary functions to get element mode
inline auto elementModeGl(const SurfaceDisplaySettings&) {
    return GL_TRIANGLES;
}
inline auto elementModeGl(const LineDisplaySettings&) {
    return GL_LINES;
}

inline auto displayGeometryType(const MembraneDisplaySettings&) {
    return DisplayGeometryType::surface;
}
inline auto displayGeometryType(const FilamentDisplaySettings& settings) {
    return settings.pathMode == FilamentDisplaySettings::PathMode::line ?
        DisplayGeometryType::line :
        DisplayGeometryType::surface;
}
inline auto displayGeometryType(const LinkerDisplaySettings& settings) {
    return settings.pathMode == LinkerDisplaySettings::PathMode::line ?
        DisplayGeometryType::line :
        DisplayGeometryType::surface;
}
inline auto displayGeometryType(const AuxLineDisplaySettings&) {
    return DisplayGeometryType::line;
}




} // namespace medyan::visual

#endif
