#ifndef MEDYAN_Visual_Render_Sphere_Hpp
#define MEDYAN_Visual_Render_Sphere_Hpp

#include <array>
#include <cstdint> // uint_fast8_t
#include <tuple>
#include <vector>

#include "Util/Math/Vec.hpp"

namespace visual {

// This function transforms a bead to GL_TRIANGLES compatible vertex coord list
// and index list.
template< typename Float = double >
struct SphereUv {
    using CoordType = mathfunc::Vec< 3, Float >;

    Float radius; // of circumscribing sphere
    std::uint_fast8_t longitudeSegs; // at least 3, about twice of latitudeSegs
    std::uint_fast8_t latitudeSegs; // at least 2

    // This function transforms a point in space to the mesh of a ball
    // Parameters
    //   - coords:  the container where the coordinates of beads can be found
    //   - indices: the bead indices on the path in the coords container
    auto generate(const CoordType& coord) const {
        using namespace mathfunc;
        using std::size_t;

        constexpr CoordType e0 = { 1.0, 0.0, 0.0 };
        constexpr CoordType e1 = { 0.0, 1.0, 0.0 };
        const CoordType e2 = cross(e0, e1);
        // The poles are at  coord +/- e2 * radius

        const size_t numVertices = longitudeSegs * (latitudeSegs - 1) + 2; // 2 means polar vertices
        const size_t numTriangles = longitudeSegs * (2 * latitudeSegs - 2);
        const auto dPhi = 2 * M_PI / longitudeSegs;
        const auto dTheta = M_PI / latitudeSegs;

        // Results
        std::vector< CoordType > vertices(numVertices); // poles at index nv - 2 and nv - 1
        std::vector< std::array< size_t, 3 > > triInd(numTriangles);

        // Calculate vertices
        for(std::uint_fast8_t nTheta = 1; nTheta < latitudeSegs; ++nTheta) {
            for(std::uint_fast8_t nPhi = 0; nPhi < longitudeSegs; ++nPhi) {
                vertices[(nTheta - 1) * longitudeSegs + nPhi]
                    = coord + radius * (
                        e0 * static_cast<Float>(std::sin(nTheta * dTheta) * std::cos(nPhi * dPhi)) +
                        e1 * static_cast<Float>(std::sin(nTheta * dTheta) * std::sin(nPhi * dPhi)) +
                        e2 * static_cast<Float>(std::cos(nTheta * dTheta))
                    );
            }
        }
        vertices[numVertices - 2] = coord + radius * e2; // Nouth pole
        vertices[numVertices - 1] = coord - radius * e2; // Sorth pole

        // Register triangles (GL_TRIANGLES)
        // from North pole to South pole
        for(std::uint_fast8_t nTheta = 0; nTheta < latitudeSegs; ++nTheta) {
            const size_t numPrev = (nTheta == 0 ? 0 : (2 * nTheta - 1) * longitudeSegs);
            for(std::uint_fast8_t nPhi = 0; nPhi < longitudeSegs; ++nPhi) {
                if(nTheta == 0) {
                    triInd[numPrev + nPhi][0] = numVertices - 2; // North pole
                    triInd[numPrev + nPhi][1] = nTheta * longitudeSegs + nPhi;
                    triInd[numPrev + nPhi][2] = nTheta * longitudeSegs + (nPhi + 1) % longitudeSegs;
                }
                else if(nTheta == latitudeSegs - 1) {
                    triInd[numPrev + nPhi][0] = (nTheta - 1) * longitudeSegs + nPhi;
                    triInd[numPrev + nPhi][1] = numVertices - 1; // South pole
                    triInd[numPrev + nPhi][2] = (nTheta - 1) * longitudeSegs + (nPhi + 1) % longitudeSegs;
                }
                else {
                    triInd[numPrev + 2 * nPhi][0]     = (nTheta - 1) * longitudeSegs + nPhi;
                    triInd[numPrev + 2 * nPhi][1]     =  nTheta      * longitudeSegs + nPhi;
                    triInd[numPrev + 2 * nPhi][2]     = (nTheta - 1) * longitudeSegs + (nPhi + 1) % longitudeSegs;
                    triInd[numPrev + 2 * nPhi + 1][0] = (nTheta - 1) * longitudeSegs + (nPhi + 1) % longitudeSegs;
                    triInd[numPrev + 2 * nPhi + 1][1] =  nTheta      * longitudeSegs + nPhi;
                    triInd[numPrev + 2 * nPhi + 1][2] =  nTheta      * longitudeSegs + (nPhi + 1) % longitudeSegs;
                }
            }
        }

        return std::make_tuple(vertices, triInd);
    }

};

} // namespace visual

#endif
