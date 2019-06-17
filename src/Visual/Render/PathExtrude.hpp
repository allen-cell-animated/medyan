#ifndef MEDYAN_Visual_Render_PathExtrude_Hpp
#define MEDYAN_Visual_Render_PathExtrude_Hpp

#include <cstdint> // uint_fast8_t
#include <tuple>
#include <vector>

#include "Util/Math/Vec.hpp"

namespace visual {

// This function transforms a path to GL_TRIANGLES compatible vertex coord list
// and index list.
template< typename Float = double >
struct PathExtrude {
    using CoordType = mathfunc::Vec< 3, Float >;

    Float radius; // of circumscribing circle
    std::uint_fast8_t sides;

    // This function transforms a path to a mesh of tubes
    // Parameters
    //   - coords:  the container where the coordinates of beads can be found
    //   - indices: the bead indices on the path in the coords container
    template< typename CoordContainer, typename IndexContainer >
    auto generate(const CoordContainer& coords, const IndexContainer& indices) const {
        using namespace mathfunc;
        using std::size_t;

        constexpr CoordType a0 { 1.0, 0.0, 0.0 };
        constexpr CoordType a1 { 0.0, 1.0, 0.0 }; // Unused unless first segment is parallel to a0

        const size_t numVertices = indices.size();
        const size_t numTubeVertices = sides * numVertices;

        // Results
        VecArray< 3, Float > vertices;
        std::vector< unsigned > triInd;

        if(indices.size() < 2)
            return std::make_tuple(vertices, triInd);

        // Locate first circle
        auto seg = normalizedVector(coords[indices[1]] - coords[indices[0]]);
        auto n0 = cross(seg, a0);
        if(magnitude2(n0) == static_cast< Float >(0.0)) n0 = cross(seg, a1);
        normalize(n0);
        auto n1 = normalizedVector(cross(n0, seg));
        auto segn = seg;

        for(size_t j = 0; j < sides; ++j) {
            const auto a = j * 2 * M_PI / sides;
            const auto point = coords[indices[0]] + n0 * (radius * std::cos((Float)a)) + n1 * (radius * std::sin((Float)a));
            vertices.push_back(point);
        }

        // Propagate circles
        for(size_t i = 1; i < numVertices; ++i) {
            segn = normalizedVector(i == numVertices - 1 ? coords[indices[i]] - coords[indices[i-1]] : coords[indices[i+1]] - coords[indices[i]]);
            const auto t = normalizedVector(0.5 * (seg + segn));
            for(size_t j = 0; j < sides; ++j) {
                // Solve for p_new given:
                //   dot(p_new - coords[indices[i]], t) = 0
                //   p_new = x * seg + pp
                const auto pp = vertices[(i-1) * sides + j];
                const auto x = dot(coords[indices[i]] - pp, t) / dot(seg, t);
                vertices.push_back(x * seg + pp);
            }

            seg = segn;
        }

        // Make indices (GL_TRIANGLES)
        triInd.resize((numVertices - 1) * (6 * sides));
        for(size_t i = 0; i < numVertices - 1; ++i) {
            for(size_t j = 0; j < sides; ++j) {
                // First triangle
                triInd[i * (6 * sides) + 6 * j    ] = (i    ) * sides + j;
                triInd[i * (6 * sides) + 6 * j + 1] = (i + 1) * sides + j;
                triInd[i * (6 * sides) + 6 * j + 2] = (i    ) * sides + (j + 1) % sides;
                // Second triangle
                triInd[i * (6 * sides) + 6 * j + 3] = (i    ) * sides + (j + 1) % sides;
                triInd[i * (6 * sides) + 6 * j + 4] = (i + 1) * sides + j;
                triInd[i * (6 * sides) + 6 * j + 5] = (i + 1) * sides + (j + 1) % sides;
            }
        }

        return std::make_tuple(vertices, triInd);
    }

};

} // namespace visual

#endif
