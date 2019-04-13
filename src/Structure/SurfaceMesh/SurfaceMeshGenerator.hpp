#ifndef MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp
#define MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp

#include <algorithm> // max, min
#include <array>
#include <cstdint> // uint_fast8_t
#include <type_traits> // integral_constant, enable_if
#include <vector>

#include "util/math/vec.hpp"

namespace mesh_gen {

namespace indexer {

using uif8 = std::uint_fast8_t;

//-------------------------------------------------------------------------
// Indexing in a tetrahedron:
//
// The vertices of the tetrahedra are v0, v1, v2, v3, which must satisfy
//     r01 x r12 . r23 > 0 (specifically, the cuboid volume)
//
// The edges in a tetrahedra are in the following order:
//     01, 02, 03, 12, 13, 23
//
// The faces in a tetrahedra are ordered the same as the opposing vertex
//-------------------------------------------------------------------------

namespace internal {

    template< uif8 v0, uif8 v1 > struct TetraEdgeIndexFromVertexIndex;
    template<> struct TetraEdgeIndexFromVertexIndex< 0, 1 > : std::integral_constant< uif8, 0 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 0, 2 > : std::integral_constant< uif8, 1 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 0, 3 > : std::integral_constant< uif8, 2 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 1, 2 > : std::integral_constant< uif8, 3 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 1, 3 > : std::integral_constant< uif8, 4 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 2, 3 > : std::integral_constant< uif8, 5 > {};
    template< uif8 v0, uif8 v1 >
    constexpr uif8 tetraEdgeIndexFromVertexIndex = TetraEdgeIndexFromVertexIndex< v0, v1 >::value;

    // 4 is invalid value
    constexpr uif8 tetraFaceIndexFromEdgeIndexTable[6][6] {
        4, 3, 2, 3, 2, 4,
        3, 4, 1, 3, 4, 1,
        2, 1, 4, 4, 2, 1,
        3, 3, 4, 4, 0, 0,
        2, 4, 2, 0, 4, 0,
        4, 1, 1, 0, 0, 4
    };
    template< uif8 e0, uif8 e1, std::enable_if_t< e0 != e1 && e0 + e1 != 5 >* = nullptr >
    constexpr uif8 getTetraFaceIndexFromEdgeIndex() {
        return tetraFaceIndexFromEdgeIndexTable[e0][e1];
    }
    template< uif8 e0, uif8 e1 >
    constexpr uif8 tetraFaceIndexFromEdgeIndex = getTetraFaceIndexFromEdgeIndex< e0, e1 >();

} // namespace internal

struct TetraTriangleIntersectionIndex {
    uif8 vertexIndex[3][2];
    uif8 edgeIndex[3];
    uif8 faceIndex[3]; // (e2-e0, e0-e1, e1-e2)
};
template< uif8 v00, uif8 v01, uif8 v10, uif8 v11, uif8 v20, uif8 v21 >
constexpr auto getTetraTriangleIntersectionIndex() {
    constexpr uif8 e0 = internal::tetraEdgeIndexFromVertexIndex< v00, v01 >;
    constexpr uif8 e1 = internal::tetraEdgeIndexFromVertexIndex< v10, v11 >;
    constexpr uif8 e2 = internal::tetraEdgeIndexFromVertexIndex< v20, v21 >;
    return TetraTriangleIntersectionIndex {
        { {v00, v01}, {v10, v11}, {v20, v21} },
        { e0, e1, e2 },
        {
            internal::tetraFaceIndexFromEdgeIndex< e2, e0 >,
            internal::tetraFaceIndexFromEdgeIndex< e0, e1 >,
            internal::tetraFaceIndexFromEdgeIndex< e1, e2 >
        }
    };
}

struct TetraQuadIntersectionIndex {
    uif8 vertexIndex[4][2];
    uif8 edgeIndex[4];
    uif8 faceIndex[4]; // (e3-e0, e0-e1, e1-e2, e2-e3)
};
template< uif8 v00, uif8 v01, uif8 v10, uif8 v11, uif8 v20, uif8 v21, uif8 v30, uif8 v31 >
constexpr auto getTetraQuadIntersectionIndex() {
    constexpr uif8 e0 = internal::tetraEdgeIndexFromVertexIndex< v00, v01 >;
    constexpr uif8 e1 = internal::tetraEdgeIndexFromVertexIndex< v10, v11 >;
    constexpr uif8 e2 = internal::tetraEdgeIndexFromVertexIndex< v20, v21 >;
    constexpr uif8 e3 = internal::tetraEdgeIndexFromVertexIndex< v30, v31 >;
    return TetraQuadIntersectionIndex {
        { {v00, v01}, {v10, v11}, {v20, v21}, {v30, v31} },
        { e0, e1, e2, e3 },
        {
            internal::tetraFaceIndexFromEdgeIndex< e3, e0 >,
            internal::tetraFaceIndexFromEdgeIndex< e0, e1 >,
            internal::tetraFaceIndexFromEdgeIndex< e1, e2 >,
            internal::tetraFaceIndexFromEdgeIndex< e2, e3 >
        }
    };
}

} // namespace indexer

template< typename Float = double >
class MarchingTetrahedraGenerator {
public:
    using small_size_t = std::uint_fast8_t;

    struct TetraEdgeData {
        Float position; // range [0.0, 1.0), might be adjusted by minimum shift from either end
        std::size_t vertexIdxInMesh;
        bool hasIntersection = false;
    };
    struct TetraEdgeIndex {
        std::size_t edgeIndex; // index in the cuboid system
        bool flipped;          // whether edge has opposite direction specified in tetra
    };

    struct TetraFaceData {
        std::size_t firstHalfEdgeIdxInMesh;
        bool hasHalfEdge = false;
    };

    MarchingTetrahedraGenerator(
        Float cubeSize,
        const mathfunc::Vec< 3, Float >& boundingBoxOrigin,
        const std::array< std::size_t, 3 >& numCubes
    ) : _cuboidSize{ cubeSize, cubeSize, cubeSize },
        _numCuboids(numCubes),
        _boundingBoxOrigin(boundingBoxOrigin),
        _distValue(_getVertexListSize()),
        _edgeData(_getEdgeListSize())
    { }

    // The function that generates the meshwork
    // Template parameters:
    //   - Mesh: the type of the surface meshwork (with attributes)
    //   - FieldFunc: a functor that takes a point and returns given how far it is outside the surface
    template< typename Mesh, typename FieldFunc >
    void operator()(Mesh& mesh, FieldFunc&& func) {
        using std::size_t;

        // Precompute phase
        for(size_t nx = 0; nx <= _numCuboids[0]; ++nx)
            for(size_t ny = 0; ny <= _numCuboids[1]; ++ny)
                for(size_t nz = 0; nz <= _numCuboids[2]; ++nz)
                    _distValue[_getVertexIdx(nx, ny, nz)] = func(_vertexCoordinate(nx, ny, nz));

        // Generating triangles
        for(size_t nx = 0; nx < _numCuboids[0]; ++nx)
            for(size_t ny = 0; ny < _numCuboids[1]; ++ny)
                for(size_t nz = 0; nz < _numCuboids[2]; ++nz)
                    for(small_size_t tetIdx = 0; tetIdx < _numTetrahedraPerCuboid; ++tetIdx)
                        calcTetra(mesh, nx, ny, nz, tetIdx);

    } // void operator()(...)

    // Helper function: computing intersection for a tetrahedron
    template< typename Mesh >
    void calcTetra(Mesh& mesh, std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx) {
        using namespace std;
        using namespace mathfunc;

        // switch sign of 4 switches and decide which procedure to follow (0.0 count as positive)
        // all pos / neg: nothing
        // one pos: gen triangle pointing pos
        // one neg: gen triangle pointing neg
        // 2 pos 2 neg: gen 2 triangles

        const array< size_t, 4 > vertexIndices {
            _getVertexIdxInTetra(nx, ny, nz, tetIdx, 0),
            _getVertexIdxInTetra(nx, ny, nz, tetIdx, 1),
            _getVertexIdxInTetra(nx, ny, nz, tetIdx, 2),
            _getVertexIdxInTetra(nx, ny, nz, tetIdx, 3)
        };
        const array< Float, 4 > vertexValues {
            _distValue[vertexIndices[0]],
            _distValue[vertexIndices[1]],
            _distValue[vertexIndices[2]],
            _distValue[vertexIndices[3]]
        };
        const array< Vec< 3, Float >, 4 > vertexCoords {
            _vertexCoordinateIdxInTetra(nx, ny, nz, tetIdx, 0),
            _vertexCoordinateIdxInTetra(nx, ny, nz, tetIdx, 1),
            _vertexCoordinateIdxInTetra(nx, ny, nz, tetIdx, 2),
            _vertexCoordinateIdxInTetra(nx, ny, nz, tetIdx, 3)
        };
        const array< TetraEdgeIndex, 6 > edgeIndices {
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 0),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 1),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 2),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 3),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 4),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 5)
        };
        const array< size_t, 4 > faceIndices {
            _getFaceIdxInTetra(nx, ny, nz, tetIdx, 0),
            _getFaceIdxInTetra(nx, ny, nz, tetIdx, 1),
            _getFaceIdxInTetra(nx, ny, nz, tetIdx, 2),
            _getFaceIdxInTetra(nx, ny, nz, tetIdx, 3)
        };

        const auto easyIntersect = [&](small_size_t eIdx, small_size_t v0Idx, small_size_t v1Idx) {
            newIntersect(
                mesh,
                edgeIndices[eIdx],
                vertexValues[v0Idx], vertexValues[v1Idx],
                vertexCoords[v0Idx], vertexCoords[v1Idx]
            );
        };
        const auto easyTriangle = [&](indexer::TetraTriangleIntersectionIndex idx) {
            // Create intersections
            easyIntersect(idx.edgeIndex[0], idx.vertexIndex[0][0], idx.vertexIndex[0][1]);
            easyIntersect(idx.edgeIndex[1], idx.vertexIndex[1][0], idx.vertexIndex[1][1]);
            easyIntersect(idx.edgeIndex[2], idx.vertexIndex[2][0], idx.vertexIndex[2][1]);

            // Create a new triangle
            sewTriangle(
                mesh,
                edgeIndices[idx.edgeIndex[0]], edgeIndices[idx.edgeIndex[1]], edgeIndices[idx.edgeIndex[2]],
                faceIndices[idx.faceIndex[0]], faceIndices[idx.faceIndex[1]], faceIndices[idx.faceIndex[2]]
            );
        };
        const auto easyQuad = [&](indexer::TetraQuadIntersectionIndex idx) {
            // Create intersections
            easyIntersect(idx.edgeIndex[0], idx.vertexIndex[0][0], idx.vertexIndex[0][1]);
            easyIntersect(idx.edgeIndex[1], idx.vertexIndex[1][0], idx.vertexIndex[1][1]);
            easyIntersect(idx.edgeIndex[2], idx.vertexIndex[2][0], idx.vertexIndex[2][1]);
            easyIntersect(idx.edgeIndex[3], idx.vertexIndex[3][0], idx.vertexIndex[3][1]);

            // Create a new quad
            sewQuad(
                mesh,
                edgeIndices[idx.edgeIndex[0]], edgeIndices[idx.edgeIndex[1]],
                edgeIndices[idx.edgeIndex[2]], edgeIndices[idx.edgeIndex[3]],
                faceIndices[idx.faceIndex[0]], faceIndices[idx.faceIndex[1]],
                faceIndices[idx.faceIndex[2]], faceIndices[idx.faceIndex[3]]
            );
        };

        small_size_t cond = 0; // 0b x x x x <-- each bit is 1 if sign is pos, 0 if sign is neg
        for(small_size_t i = 0; i < 4; ++i) {
            cond << 1;
            cond |= (vertexValues[i] >= 0.0 ? 1 : 0);
        }

        switch(cond) {
        case 0b0000:
        case 0b1111:
            break;

        case 0b0001:
            // intersects 03, 13, 23, triangle (03, 13, 23)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                0, 3, 1, 3, 2, 3
            >());
            break;
        case 0b1110:
            // intersects 23, 13, 03, triangle (23, 13, 03)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                2, 3, 1, 3, 0, 3
            >());
            break;

        case 0b0010:
            // intersects 12, 02, 23, triangle (12, 02, 23)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                1, 2, 0, 2, 2, 3
            >());
            break;
        case 0b1101:
            // intersects 23, 02, 12, triangle (23, 02, 12)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                2, 3, 0, 2, 1, 2
            >());
            break;

        case 0b0100:
            // intersects 12, 13, 01, triangle (12, 13, 01)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                1, 2, 1, 3, 0, 1
            >());
            break;
        case 0b1011:
            // intersects 01, 13, 12, triangle (01, 13, 12)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                0, 1, 1, 3, 1, 2
            >());
            break;

        case 0b1000:
            // intersects 03, 02, 01, triangle (03, 02, 01)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                0, 3, 0, 2, 0, 1
            >());
            break;
        case 0b0111:
            // intersects 01, 02, 03, triangle (01, 02, 03)
            easyTriangle(indexer::getTetraTriangleIntersectionIndex<
                0, 1, 0, 2, 0, 3
            >());
            break;

        case 0b0011:
            // intersects 02, 03, 13, 12, 2 triangles
            easyQuad(indexer::getTetraQuadIntersectionIndex<
                0, 2, 0, 3, 1, 3, 1, 2
            >());
            break;
        case 0b1100:
            // intersects 12, 13, 03, 02, 2 triangles
            easyQuad(indexer::getTetraQuadIntersectionIndex<
                1, 2, 1, 3, 0, 3, 0, 2
            >());
            break;

        case 0b0101:
            // intersects 01, 12, 23, 03, 2 triangles
            easyQuad(indexer::getTetraQuadIntersectionIndex<
                0, 1, 1, 2, 2, 3, 0, 3
            >());
            break;
        case 0b1010:
            // intersects 03, 23, 12, 01, 2 triangles
            easyQuad(indexer::getTetraQuadIntersectionIndex<
                0, 3, 2, 3, 1, 2, 0, 1
            >());
            break;

        case 0b0110:
            // intersects 01, 02, 23, 13, 2 triangles
            easyQuad(indexer::getTetraQuadIntersectionIndex<
                0, 1, 0, 2, 2, 3, 1, 3
            >());
            break;
        case 0b1001:
            // intersects 13, 23, 02, 01, 2 triangles
            easyQuad(indexer::getTetraQuadIntersectionIndex<
                1, 3, 2, 3, 0, 2, 0, 1
            >());
            break;

        } // switch(cond)
    } // void calcTetra(...)

    // Helper function: add a new intersection on edge
    // value0 and value1 should be the sequence in tetrahedra. This function will manage the flipping.
    // The caller must ensure that value0 and value1 are in opposite signs
    template< typename Mesh >
    void newIntersect(
        Mesh& mesh,
        const TetraEdgeIndex& idx,
        Float value0, Float value1,
        const mathfunc::Vec<3, Float>& coord0, const mathfunc::Vec<3, Float>& coord1
    ) {
        struct InsertionCoordinate {
            mathfunc::Vec< 3, Float > coord;
            auto coordinate(Mesh& mesh, size_t vi) const { return coord; }
        };

        if(!_edgeData[idx.edgeIndex].hasIntersection) {
            // Compute position
            Float pos = value0 / (value0 - value1);
            pos = std::min(1 - _minPositionShift, std::max(_minPositionShift, pos));
            auto newCoord = coord0 * pos + coord1 * (1 - pos);

            // Add new vertex
            const auto vi = typename Mesh::NewVertexInsertion{} (mesh, InsertionCoordinate{ newCoord });

            // Update edge data
            _edgeData[idx.edgeIndex].position = (idx.flipped ? pos : 1 - pos);
            _edgeData[idx.edgeIndex].hasIntersection = true;
            _edgeData[idx.edgeIndex].vertexIdxInMesh = vi;
        } // else do nothing
    }

    // Helper function: sewing a new triangle to the surface
    // The tetra edges specified should be arranged so that the direction points outward.
    template< typename Mesh >
    void sewTriangle(
        Mesh& mesh,
        small_size_t e0Idx, small_size_t e1Idx, small_size_t e2Idx,
        small_size_t f0Idx, small_size_t f1Idx, small_size_t f2Idx
    ) const {
        using std::size_t;

        const auto newHalfEdgeIndices = typename Mesh::NewTrianglePatch{} (
            mesh,
            {_edgeData[e0Idx].vertexInMesh, _edgeData[e1Idx].vertexInMesh, _edgeData[e2Idx].vertexInMesh},
            {_faceData[f0Idx].hasHalfEdge, _faceData[f1Idx].hasHalfEdge, _faceData[f2Idx].hasHalfEdge},
            {_faceData[f0Idx].firstHalfEdgeIdxInMesh, _faceData[f1Idx].firstHalfEdgeIdxInMesh, _faceData[f2Idx].firstHalfEdgeIdxInMesh}
        );

        _faceData[f0Idx].hasHalfEdge = true;
        _faceData[f1Idx].hasHalfEdge = true;
        _faceData[f2Idx].hasHalfEdge = true;
        _faceData[f0Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices[0];
        _faceData[f1Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices[1];
        _faceData[f2Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices[2];
    }

    // Helper function: sewing a (skew) quadrilateral to the surface
    // The tetra edges specified should be arranged so that the direction points outward.
    template< typename Mesh >
    void sewQuad(
        Mesh& mesh,
        small_size_t e0Idx, small_size_t e1Idx, small_size_t e2Idx, small_size_t e3Idx,
        small_size_t f0Idx, small_size_t f1Idx, small_size_t f2Idx, small_size_t f3Idx
    ) const {
        using std::size_t;

        // Add triangle with vertices on edges 012 and edges 230
        //
        //     e3 ------- f0 ------- e0
        //     |                   / |
        //     |                /    |
        //     |             /       |
        //     f3         /          f1
        //     |       /             |
        //     |    /                |
        //     | /                   |
        //     e2 ------- f2 ------- e1
        const auto newHalfEdgeIndices0 = typename Mesh::NewTrianglePatch{} (
            mesh,
            {_edgeData[e0Idx].vertexInMesh, _edgeData[e1Idx].vertexInMesh, _edgeData[e2Idx].vertexInMesh},
            {false, _faceData[f1Idx].hasHalfEdge, _faceData[f2Idx].hasHalfEdge},
            {0, _faceData[f1Idx].firstHalfEdgeIdxInMesh, _faceData[f2Idx].firstHalfEdgeIdxInMesh}
        );
        const auto newHalfEdgeIndices1 = typename Mesh::NewTrianglePatch{} (
            mesh,
            {_edgeData[e2Idx].vertexInMesh, _edgeData[e3Idx].vertexInMesh, _edgeData[e0Idx].vertexInMesh},
            {true, _faceData[f3Idx].hasHalfEdge, _faceData[f0Idx].hasHalfEdge},
            {newHalfEdgeIndices0[0], _faceData[f3Idx].firstHalfEdgeIdxInMesh, _faceData[f0Idx].firstHalfEdgeIdxInMesh}
        );

        _faceData[f0Idx].hasHalfEdge = true;
        _faceData[f1Idx].hasHalfEdge = true;
        _faceData[f2Idx].hasHalfEdge = true;
        _faceData[f3Idx].hasHalfEdge = true;
        _faceData[f0Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices1[2];
        _faceData[f1Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices0[1];
        _faceData[f2Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices0[2];
        _faceData[f3Idx].firstHalfEdgeIdxInMesh = newHalfEdgeIndices1[1];
    }

private:
    // Constants
    //-------------------------------------------------------------------------
    // Indexing in cuboid:
    //
    // Number of vertices: (nx + 1)(ny + 1)(nz + 1)
    //
    // Tetrahedra in a cuboid is ordered in the following order (labeled using r's):
    //     ijk, (i+j)(-i)(k+i), jki, (j+k)(-j)(i+j), kij, (k+i)(-k)(j+k)
    //
    // Edges in a cuboid is simply indexed by 0bxyz(edge direction) - 1
    //
    // Faces in a cuboid is ordered for each tetrahedron, with faces opposing
    // vertices
    //     v1 (inside cuboid), v3 (on the surface)
    // The opposing face for v0 never belongs to this cuboid, while the opposing
    // face for v2 is the same as the opposing face of v1 in the previous tetra.
    //-------------------------------------------------------------------------
    static constexpr small_size_t _numTetrahedraPerCuboid = 6;
    static constexpr small_size_t _numTetraFacesPerCuboid = 12; // 6 on 3 faces and 6 inside the cuboid
    static constexpr small_size_t _numTetraEdgesPerCuboid = 7; // index = 0bxyz - 1

    // Local vertex index [local tetra idx (6)][vertex idx (4)]
    static constexpr small_size_t _tetraVertexLocalIndex[6][4] {
        0b000, 0b100, 0b110, 0b111,
        0b000, 0b110, 0b010, 0b111,
        0b000, 0b010, 0b011, 0b111,
        0b000, 0b011, 0b001, 0b111,
        0b000, 0b001, 0b101, 0b111,
        0b000, 0b101, 0b100, 0b111
    };
    // Local edge index [local tetra idx (6)][edge idx (6)]
    // Value:  0b  xyz                  d                      xyz
    //             ^ starting position  ^ same(0)/diff(1) dir  ^ edge direction (positive)
    // Real index in cuboid = 0bxyz(direction) - 1
    static constexpr small_size_t _tetraEdgeLocalIndex[6][6] {
        0b0000100, 0b0000110, 0b0000111, 0b1000010, 0b1000011, 0b1100001,
        0b0000110, 0b0000010, 0b0000111, 0b1101100, 0b1100001, 0b0100101,
        0b0000010, 0b0000011, 0b0000111, 0b0100001, 0b0100101, 0b0110100,
        0b0000011, 0b0000001, 0b0000111, 0b0111010, 0b0110100, 0b0010110,
        0b0000001, 0b0000101, 0b0000111, 0b0010100, 0b0010110, 0b1010010,
        0b0000101, 0b0000100, 0b0000111, 0b1011001, 0b1010010, 0b1000011
    };
    // Local face index [local tetra idx (6)][face idx (4)]
    // Value:    0b  xyz                aaaa
    //               ^ cuboid position  ^ face id (0000(0) - 1011(11))
    static constexpr small_size_t _tetraFaceLocalIndex[6][4] {
        0b1000101, 0b0000000, 0b0001010, 0b0000001,
        0b0101011, 0b0000010, 0b0000000, 0b0000011,
        0b0101001, 0b0000100, 0b0000010, 0b0000101,
        0b0010011, 0b0000110, 0b0000100, 0b0000111,
        0b0010001, 0b0001000, 0b0000110, 0b0001001,
        0b1000111, 0b0001010, 0b0001000, 0b0001011
    };

    // Parameters
    Float                           _minPositionShift = 1e-2; // The position will have a minimal shift from either end
    mathfunc::Vec< 3, Float >       _cuboidSize;
    std::array< std::size_t, 3 >    _numCuboids;
    mathfunc::Vec< 3, Float >       _boundingBoxOrigin;

    // Internal data
    std::vector< Float > _distValue; // using vertex indexing
    std::vector< TetraEdgeData > _edgeData; // using edge indexing
    std::vector< TetraFaceData > _faceData; // using face indexing

    // Indexers
    auto _getVertexListSize() const {
        return (_numCuboids[0] + 1) * (_numCuboids[1] + 1) * (_numCuboids[2] + 1);
    }
    auto _getVertexIdx(std::size_t nx, std::size_t ny, std::size_t nz) const {
        std::size_t res = nx;
        res *= _numCuboids[0] + 1; res += ny;
        res *= _numCuboids[1] + 1; res += nz;
        return res;
    }
    auto _getVertexIdxInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t vtxIdx) const {
        const small_size_t i = _tetraVertexLocalIndex[tetIdx][vtxIdx];
        return _getVertexIdx(
            nx + (i >> 2) & 1,
            ny + (i >> 1) & 1,
            nz +  i       & 1
        );
    }
    auto _getCubeListSize() const {
        return _numCuboids[0] * _numCuboids[1] * _numCuboids[2];
    }
    auto _getCubeIdx(std::size_t nx, std::size_t ny, std::size_t nz) const {
        std::size_t res = nx;
        res *= _numCuboids[0]; res += ny;
        res *= _numCuboids[1]; res += nz;
        return res;
    }
    auto _getEdgeListSize() const {
        return _numTetraEdgesPerCuboid * _getCubeListSize();
    }
    auto _getEdgeIdx(std::size_t nx, std::size_t ny, std::size_t nz, std::size_t edgeIdx) const {
        return _numTetraEdgesPerCuboid * _getCubeIdx(nx, ny, nz) + edgeIdx;
    }
    auto _getEdgeIdxInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t edgeIdx) const {
        const small_size_t i = _tetraEdgeLocalIndex[tetIdx][edgeIdx];
        const small_size_t idxInCuboid = i & 0b111 - 1;
        const bool flipped = (i >> 3) & 1;
        return TetraEdgeIndex {
            _getEdgeIdx(
                nx + (i >> 6) & 1,
                ny + (i >> 5) & 1,
                nz + (i >> 4) & 1,
                idxInCuboid
            ),
            flipped
        };
    }
    auto _getFaceIdx(std::size_t nx, std::size_t ny, std::size_t nz, std::size_t faceIdx) const {
        return _numTetraFacesPerCuboid * _getCubeIdx(nx, ny, nz) + faceIdx;
    }
    auto _getFaceIdxInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t faceIdx) const {
        const small_size_t i = _tetraFaceLocalIndex[tetIdx][faceIdx];
        const small_size_t idxInCuboid = i & 0b1111;
        return _getFaceIdx(
            nx + (i >> 6) & 1,
            ny + (i >> 5) & 1,
            nz + (i >> 4) & 1,
            idxInCuboid
        );
    }

    // Coordinates
    auto _vertexCoordinate(std::size_t nx, std::size_t ny, std::size_t nz) const {
        return mathfunc::Vec< 3, Float > {
            _boundingBoxOrigin[0] + _cuboidSize[0] * nx,
            _boundingBoxOrigin[1] + _cuboidSize[1] * ny,
            _boundingBoxOrigin[2] + _cuboidSize[2] * nz
        };
    }
    auto _vertexCoordinateIdxInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t vtxIdx) const {
        const small_size_t i = _tetraVertexLocalIndex[tetIdx][vtxIdx];
        return _vertexCoordinate(
            nx + (i >> 2) & 1,
            ny + (i >> 1) & 1,
            nz +  i       & 1
        );
    }
};

// static variable definition (declaration in C++17, and should be removed)
template< typename Float >
constexpr typename MarchingTetrahedraGenerator<Float>::small_size_t
MarchingTetrahedraGenerator<Float>::_tetraVertexLocalIndex[6][4];
template< typename Float >
constexpr typename MarchingTetrahedraGenerator<Float>::small_size_t
MarchingTetrahedraGenerator<Float>::_tetraEdgeLocalIndex[6][6];
template< typename Float >
constexpr typename MarchingTetrahedraGenerator<Float>::small_size_t
MarchingTetrahedraGenerator<Float>::_tetraFaceLocalIndex[6][4];

} // namespace mesh_gen

#endif
