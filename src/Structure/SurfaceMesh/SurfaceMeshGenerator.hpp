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

    struct TetraFaceData {
        std::size_t firstHalfEdgeIdxInMesh;
        bool hasHalfEdge;
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
                    ;
    }

    // Helper function: computing intersection for a tetrahedron
    void calcTetra(std::size_t nx ny nz, small_size_t tetIdx) const {
        // switch sign of 4 switches and decide which procedure to follow (0.0 count as positive)
        // all pos / neg: nothing
        // one pos: gen triangle pointing pos
        // one neg: gen triangle pointing neg
        // 2 pos 2 neg: gen 2 triangles
        get_vtx_idx from (nx ny nz and tetIdx);

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
    void sewTriangle(
        std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx,
        small_size_t e0Idx, small_size_t e1Idx, small_size_t e2Idx
    ) const {
        get face idx from nx ny nz and tetIdx, e Indices

        tell "mesh" to add a new triangle: (v0, v1, v2, he0op, he0opIdx, he1op, he1opIdx, he2op, he2opIdx)
        update face data =
        for i = 0:3 get face data:
            if face has an halfedge (ref target vertex, and edge) -> create opposite
            else -> create a single halfedge (ref target vertex) and an edge (ref complete); -> update halfedge (ref edge)
        create triangle (ref complete)
        update halfedges (ref triangle, next, prev)

        update ref for vertices
    }

private:
    // Constants
    static constexpr small_size_t _numTetrahedraPerCuboid = 6;
    static constexpr small_size_t _numTetraFacesPerCuboid = 12; // 6 on 3 faces and 6 inside the cuboid
    static constexpr small_size_t _numTetraEdgesPerCuboid = 7; // 001 010 100, 110 101 011, 111

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

    // Coordinates
    auto _vertexCoordinate(std::size_t nx, std::size_t ny, std::size_t nz) const {
        return mathfunc::Vec< 3, Float > {
            _boundingBoxOrigin[0] + _cuboidSize[0] * nx,
            _boundingBoxOrigin[1] + _cuboidSize[1] * ny,
            _boundingBoxOrigin[2] + _cuboidSize[2] * nz
        }
    }
};

} // namespace mesh_gen

#endif
