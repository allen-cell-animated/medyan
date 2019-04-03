#ifndef MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp
#define MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp

#include <array>
#include <vector>

#include "util/math/vec.hpp"

namespace mesh_gen {

template< typename Float = double >
class MarchingTetrahedraGenerator {
public:

    struct TetraEdgeData {
        Float position; // range [0.0, 1.0)
        std::size_t vertexIdxInMesh;
        bool hasIntersection;
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

private:
    static constexpr std::size_t _numTetrahedraPerCuboid = 6;
    static constexpr std::size_t _numTetraFacesPerCuboid = 12;
    static constexpr std::size_t _numTetraEdgesPerCuboid = 7; // 001 010 100, 110 101 011, 111

    mathfunc::Vec< 3, Float >       _cuboidSize;
    std::array< std::size_t, 3 >    _numCuboids;
    mathfunc::Vec< 3, Float >       _boundingBoxOrigin;

    std::vector< Float > _distValue; // using vertex indexing
    std::vector< TetraEdgeData > _edgeData; // using edge indexing

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
