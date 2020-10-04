#ifndef MEDYAN_Visual_Render_PathExtrude_Hpp
#define MEDYAN_Visual_Render_PathExtrude_Hpp

#include <array>
#include <cstdint> // uint_fast8_t
#include <tuple>
#include <vector>

#include "Util/Math/Vec.hpp"

namespace medyan::visual {

// This function transforms a path to GL_TRIANGLES compatible vertex coord list
// and index list.
template< typename Float = double >
struct PathExtrude {
    using CoordType = mathfunc::Vec< 3, Float >;

    Float radius; // of circumscribing circle
    int   sides;

    // This function transforms a path to a mesh of tubes
    // Parameters
    //   - coords:  the container where the coordinates of beads can be found
    //   - indices: the bead indices on the path in the coords container
    // Future improvements:
    //   - Introduce different vertex normals at conjunctions.
    template< typename CoordContainer, typename IndexContainer >
    auto generate(const CoordContainer& coords, const IndexContainer& indices) const {
        using namespace mathfunc;
        using std::size_t;

        constexpr CoordType a0 { 1.0, 0.0, 0.0 };
        constexpr CoordType a1 { 0.0, 1.0, 0.0 }; // Unused unless first segment is parallel to a0

        const size_t numVertices = indices.size();
        const size_t numTubeVertices = sides * numVertices;
        const size_t numTriangles = (numVertices - 1) * 2 * sides;

        // Results
        std::vector< CoordType > vertices;      vertices.reserve(numTubeVertices);
        std::vector< CoordType > vertexNormals; vertexNormals.reserve(numTubeVertices);
        std::vector< std::array< int, 3 > > triInd(numTriangles);

        if(indices.size() < 2)
            return std::make_tuple(vertices, vertexNormals, triInd);

        // Locate first circle
        CoordType seg ( normalizedVector(coords[indices[1]] - coords[indices[0]]) );
        auto n0 = cross(seg, a0);
        if(magnitude2(n0) == static_cast< Float >(0.0)) n0 = cross(seg, a1);
        normalize(n0);
        auto n1 = normalizedVector(cross(n0, seg));
        auto segn = seg;

        for(size_t j = 0; j < sides; ++j) {
            const Float a = j * 2 * M_PI / sides;
            const Float cosa = std::cos(a);
            const Float sina = std::sin(a);
            const auto point = static_cast< CoordType >(coords[indices[0]]) + n0 * (radius * cosa) + n1 * (radius * sina);
            const auto un    = n0 * cosa + n1 * sina;
            vertices.push_back(point);
            vertexNormals.push_back(un);
        }

        // Propagate circles
        for(size_t i = 1; i < numVertices; ++i) {
            segn = normalizedVector(i == numVertices - 1 ? coords[indices[i]] - coords[indices[i-1]] : coords[indices[i+1]] - coords[indices[i]]);
            const auto t = normalizedVector(seg + segn);
            const auto dot_seg_t = dot(seg, t);

            for(size_t j = 0; j < sides; ++j) {
                // Solve for p_new given:
                //   dot(p_new - coords[indices[i]], t) = 0
                //   p_new = x * seg + pp
                const auto pp = vertices[(i-1) * sides + j];
                const Float x = dot(coords[indices[i]] - pp, t) / dot_seg_t;
                vertices.push_back(x * seg + pp);
                vertexNormals.push_back(normalizedVector(vertices.back() - static_cast< CoordType >(coords[indices[i]])));
            }

            seg = segn;
        }

        // Make indices (GL_TRIANGLES)
        for(size_t i = 0; i < numVertices - 1; ++i) {
            for(size_t j = 0; j < sides; ++j) {
                // First triangle
                triInd[i * (2 * sides) + 2 * j][0] = (i    ) * sides + j;
                triInd[i * (2 * sides) + 2 * j][1] = (i + 1) * sides + j;
                triInd[i * (2 * sides) + 2 * j][2] = (i    ) * sides + (j + 1) % sides;
                // Second triangle
                triInd[i * (2 * sides) + 2 * j + 1][0] = (i    ) * sides + (j + 1) % sides;
                triInd[i * (2 * sides) + 2 * j + 1][1] = (i + 1) * sides + j;
                triInd[i * (2 * sides) + 2 * j + 1][2] = (i + 1) * sides + (j + 1) % sides;
            }
        }

        return std::make_tuple(vertices, vertexNormals, triInd);
    }

    // Provide an estimate on the number of triangles returned by "generate"
    // function.
    int estimateNumTriangles(int numVertices) const {
        if(numVertices < 2) return 0;

        const int numSegments = numVertices - 1;
        return numSegments * sides * 2;
    }

};

} // namespace medyan::visual

#endif
