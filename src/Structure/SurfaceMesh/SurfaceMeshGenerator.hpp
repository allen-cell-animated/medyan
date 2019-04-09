#ifndef MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp
#define MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp

#include <array>
#include <cstdint> // uint_fast8_t
#include <vector>

#include "util/math/vec.hpp"

namespace mesh_gen {

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
        std::size_t edgeIndex;
        bool flipped;
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
    void operator()(Mesh& mesh, FieldFunc&&) {
        using std::size_t;

        // Precompute phase
        for(size_t nx = 0; nx <= _numCuboids[0]; ++nx)
            for(size_t ny = 0; ny <= _numCuboids[1]; ++ny)
                for(size_t nz = 0; nz <= _numCuboids[2]; ++nz)
                    _distValue[_getVertexIdx(nx, ny, nz)] = FieldFunc(_vertexCoordinate(nx, ny, nz));

        // Generating triangles
        for(size_t nx = 0; nx < _numCuboids[0]; ++nx)
            for(size_t ny = 0; ny < _numCuboids[1]; ++ny)
                for(size_t nz = 0; nz < _numCuboids[2]; ++nz)
                    for(small_size_t tetIdx = 0; tetIdx < _numTetrahedraPerCuboid; ++tetIdx)
                        calcTetra(nx, ny, nz, tetIdx);
    }

    // Helper function: computing intersection for a tetrahedron
    void calcTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx) const {
        using namespace std;

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
        const array< TetraEdgeIndex, 6 > edgeIndices {
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 0),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 1),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 2),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 3),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 4),
            _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 5),
        };
        get face index from (nx ny nz and tetIdx);

        small_size_t cond = 0; // 0b x x x x <-- each bit is 1 if sign is pos, 0 if sign is neg
        update cond();

        switch(cond) {
        case 0b0000:
        case 0b1111:
            break;

        case 0b0001:
            // compute position on edges
            // compute coordinate and create new vertex (w/o edges, halfedges or triangles)
        }
    }

    // Helper function: sewing a new triangle to the surface
    template< typename Mesh >
    void sewTriangle(
        Mesh& mesh,
        std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx,
        small_size_t e0Idx, small_size_t e1Idx, small_size_t e2Idx
    ) const {
        using std::size_t;

        get face idx from nx ny nz and tetIdx, e Indices

        typename Mesh::NewTrianglePatch{} (
            mesh,
            {v??, v??, v??},
            {_faceData[??].hasHalfEdge, _faceData[??].hasHalfEdge, _faceData[??].hasHalfEdge},
            {_faceData[??].firstHalfEdgeIdxInMesh, _faceData[??].firstHalfEdgeIdxInMesh, _faceData[??].firstHalfEdgeIdxInMesh}
        );
        tell "mesh" to add a new triangle: (v0, v1, v2, he0op, he0opIdx, he1op, he1opIdx, he2op, he2opIdx)
    }

private:
    // Constants
    //-------------------------------------------------------------------------
    // Indexing in tetrahedra:
    //
    // The vertices of the tetrahedra are v0, v1, v2, v3, which must satisfy
    //     r01 x r12 . r23 > 0 (specifically, the cuboid volume)
    //
    // The edges in a tetrahedra are in the following order:
    //     01, 02, 03, 12, 13, 23
    // The faces in a tetrahedra are ordered the same as the opposing vertex
    //-------------------------------------------------------------------------
    // Indexing in cuboid:
    //
    // Tetrahedra in a cuboid is ordered in the following order (labeled using r's):
    //     ijk, (i+j)(-i)(k+i), jki, (j+k)(-j)(i+j), kij, (k+i)(-k)(j+k)
    //
    // Edges in a cuboid is simply indexed by 0bxyz(edge direction) - 1
    //
    // Faces in a cuboid is ordered the same with tetrahedra, with faces opposing vertices
    //     v1 (inside cuboid), v3 (on the surface)
    // The opposing face for v0 never belongs the this cuboid, while the opposing
    // face for v2 is the same as the opposing face of v1 in the previous tetra.
    //-------------------------------------------------------------------------
    static constexpr small_size_t _numTetrahedraPerCuboid = 6;
    static constexpr small_size_t _numTetraFacesPerCuboid = 12; // 6 on 3 faces and 6 inside the cuboid
    static constexpr small_size_t _numTetraEdgesPerCuboid = 7; // = 0bxyz - 1

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
