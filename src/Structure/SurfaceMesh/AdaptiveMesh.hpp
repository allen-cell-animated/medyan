/*

Adaptive mesh algorithm

Implementation based on
"An Adaptive Mesh Algorithm for Evolving Surfaces: Simulations of Drop Breakup and Coalescence"
by Vittorio Cristini, Jerzy Blawzdziewicz and Michael Loewenberg.

Performs mesh relaxation and topological transformation.

Ending criteria:
    1. Relaxation terminates.
    2. in a local region, number of nodes n is near n_0.
    3. no more node reconnection.

The following values are stored and updated
    - L computed on each vertex.
        - L_0 is generally universal, as an upper limit
        - L_1 uses local curvature information.
        - L_2 and above are not used currently.
    - rho computed on each vertex. rho = c_0 / (alpha * L)^2
        - c_0 is a geometric factor (2 / sqrt(3))
        - alpha_0 is resolution of the surface (ref: 0.2-0.3)
    - rho_avg computed on each vertex. (1-ring weighted avg of rho)
    - A_0 computed an each triangle. Uses rho_avg on vertices.

For mesh relaxation, the following values are also needed
    - l_0' computed on each vertex. Uses weighted sum of neighboring A_0.
    - l_0 computed on each edge. An average of l_0' on the two vertex.
    - normal computed on vertices.
    - v (velocity) on vertices, Uses l_0 and n.

For topological transformation, the following values are also needed
    - S_loc computed on each vertex
    - S_0 computed on a local region. Uses A_0 in the region.
    - n_0 computed on a local region. n_0 = n * S_loc / S_0. (fraction of border nodes will be counted to n).

We might need to redistribute the surface area globally or locally depending on the operation,
if we are using patched area energy computation (i.e. use sum of local area energy instead of
energy from sum of area)

The algorithm was not introduced explicity in the article, and we'll formalize it as follows

Init: Find per-element quantities for all elements.
Loop
    Global relaxation
    Update per-element quantities and local averaged quantities
    For all places that does not satisfy criteria 2, 3
        Local topological transformation
        Local relaxation
        Update affected per-element quantities and local averaged quantites
    End
Until all criteria are met

It would be easier if the implementation and the mesh representation are coupled with
the current surface meshwork system. But additional variables should be introduced and
it might not be appropriate to mess them up with the original structure.

Overlaying a new set of variables with the current implementation of the meshwork might
be a good idea as well.

*/

#ifndef MEDYAN_AdaptiveMesh_hpp
#define MEDYAN_AdaptiveMesh_hpp

#include <algorithm> // max min
#include <vector>

#include "MathFunctions.h"

#include "Structure/SurfaceMesh/AdaptiveMeshAttribute.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"


template< typename Mesh >
void calc_data_for_affected(Mesh& mesh) {
    for(size_t i = 0; i < mesh.numVertices(); ++i) {
        if(affected(i, lv1)) {
            calc_curv_and_gaussian_curv();
            const double l0 = <some-external-value>
            const double l1 = 1.0 / sqrt(std::max(2 * curv * curv - gaussianCurv, 0.0));
            const double l = std::min(l0, l1);
            mesh.getVertexAttribute(i).adapt.density0 = c0 * alphaInvSqr / (l*l);
        }
    }
    for(size_t i = 0; i < mesh.numVertices(); ++i) {
        if(affected(i, lv2)) {
            mesh.getVertexAttribute(i).adapt.densityAvg = 0.5 * (
                mesh.getVertexAttribute(i).adapt.density0
                + sumNeighborDensity0
            );
        }
    }

    for(size_t i = 0; i < mesh.)
    double L0 = <some-external-value>
}

template< typename Mesh > void calc_l0_all(Mesh& mesh) {
    static const double sqrt2C0 = std::sqrt(2 * c0);

    for(size_t i = 0; i < mesh.numVertices(); ++i) {
        calc_curv_and_gaussian_curv();
        const auto l0 = <some-external-value>
        const auto l1 = 1.0 / sqrt(std::max(2 * curv * curv - gaussianCurv, 0.0));
        const auto l = std::min(l0, l1);
        mesh.getVertexAttribute(i).adapt.density0 = c0 * alphaInvSqr / (l*l);
    }
    for(size_t i = 0; i < mesh.numVertices(); ++i) {
        double sumNeighborDensity0 = 0.0;
        size_t numNeighbors = 0;
        mesh.forEachHalfEdgeTargetingVertex(i, [&mesh, &sumNeighborDensity0, &numNeighbors](size_t hei) {
            sumNeighborDensity0 += mesh.getVertexAttribute(mesh.target(mesh.opposite(hei))).adapt.density0;
            ++numNeighbors;
        });
        mesh.getVertexAttribute(i).adapt.densityAvg = 0.5 * (
            mesh.getVertexAttribute(i).adapt.density0
            + sumNeighborDensity0 / numNeighbors;
        );
    }
    for(size_t i = 0; i < mesh.numTriangles(); ++i) {
        double sumDensityAvg = 0.0;
        mesh.forEachHalfEdgeInTriangle(i, [&mesh, &sumDensityAvg](size_t hei) {
            sumDensityAvg += mesh.getVertexAttribute(mesh.target(hei)).adapt.densityAvg;
        });
        mesh.getTriangleAttribute(i).adapt.area0 = 1.5 / sumDensityAvg;
    }
    for(size_t i = 0; i < mesh.numVertices(); ++i) {
        double sumNeighborWeightedArea0 = 0.0;
        size_t numNeighbors = 0;
        mesh.forEachHalfEdgeTargetingVertex(i, [&mesh, &sumNeighborWeightedArea0, &numNeighbors](size_t hei) {
            const size_t tInner = mesh.triangle(hei);
            const size_t tOuter = mesh.triangle(mesh.opposite(mesh.prev(hei)));
            sumNeighborWeightedArea0 +=
                mesh.getTriangleAttribute(tInner).adapt.area0 * 5.0 / 6.0
                + mesh.getTriangleAttribute(tOuter).adapt.area0 / 6.0;
            ++numNeighbors;
        });
        mesh.getVertexAttribute(i).adapt.l0Aux = sqrt2C0 * sumNeighborWeightedArea0 / numNeighbors;
    }
    for(size_t i = 0; i < mesh.numEdges(); ++i) {
        double suml0Aux = 0.0;
        mesh.forEachHalfEdgeInEdge(i, [&mesh, &suml0Aux](size_t hei) {
            suml0Aux += mesh.getVertexAttribute(mesh.target(hei)).adapt.l0Aux;
        });
        mesh.getEdgeAttribute(i).adapt.l0 = suml0Aux * 0.5;
    }
}

template< typename Mesh > auto compute_all_forces(const Mesh& mesh, const mathfunc::VecMut< mathfunc::Vec3 >& coords) {
    mathfunc::VecMut< mathfunc::Vec3 > forces(coords.size());
    double maxMag2F = 0.0;

    std::transform(
        mesh.getVertices().begin(), mesh.getVertices().end(),
        coords.begin(),
        forces.begin(),
        [&coords, &maxForce](const auto& v, const auto& c) {
            mathfunc::Vec3 f {};
            mesh.forEachHalfEdgeTargetingVertex(v, [&](size_t hei) {
                const auto l0 = mesh.getEdgeAttribute(mesh.edge(hei)).adapt.l0;
                const auto r = coords[mesh.target(mesh.opposite(hei))] - c;
                const auto mag = mathfunc::magnitude(r);
                f += r * (1.0 - l0 / mag);
            });

            const auto& un = v.attr.adapt.unitNormal;
            f -= un * dot(un, f); // remove normal component

            const auto mag2F = mathfunc::magnitude2(f);
            if(mag2F > maxMag2F) maxMag2F = mag2F;

            return f;
        }
    );

    return std::make_tuple(std::move(forces), maxMag2F);
}
// 2nd order Runge Kutta method on a set of coordinates.
void rk2(coord, calc_force, Float dt, Float epsilon) {
    force = calc_force(coord);
    coord_halfway = coord + force * dt / 2;
    force_halfway = calc_force(coord_halfway);
    coord += force_halfway * dt;
}

template< typename Mesh > void global_relaxation(Mesh& mesh) {
    using coordinate_type = typename Mesh::VertexAttribute::coordinate_type;

    // TODO Need epsilonSqr
    // normal is obtained from curvature computation
    double maxMag2F;
    const size_t numVertices = mesh.numVertices();

    mathfunc::VecMut< mathfunc::Vec3 > coords(numVertices);

    // Copy coordinates to a vector
    std::transform(
        mesh.getVertices().begin(), mesh.getVertices().end(),
        coords.begin(),
        [](const auto& v) { return vector2Vec<3, double>(v.attr.vertex->coordinate); }
    );
    // Compute forces
    mathfunc::VecMut< mathfunc::Vec3 > forces;

    std::tie(forces, maxMag2F) = compute_all_forces(mesh, coords);

    while(maxMag2F >= epsilonSqr) {
        // TODO: runge kutta with coords, forces

        std::tie(forces, maxMag2F) = compute_all_forces(mesh, coords);
    }
}
// general algorithm procedure
template< typename Mesh >
void adaptive_mesh(Mesh& mesh) {
    calc_l0_all();
    until(density critieria are met) {
        global_relaxation(); // Lowers the force
        while(exists density criteria violation) {
            violation = 1st violation;
            correct(violation);
            local_remesh(around corrected element);
            update_local_criteria(around corrected element);
        }

    }
}
struct MeshForce {
    const std::vector< Vertex* >& vertices;
    std::vector< char > mask;
    std::vector< size_t > activeVertexIndices;
    some_array operator()(some_array coord) {
        // The size of coord must be exactly 3 times the size of vertex array.
        for(size_t i : activeVertexIndices) {
            Vector3 p = vertices[i] -> point;
            for(Vertex* v : vertices[i] -> neighborVertices) {
                if
            }
        }
    }
};

#endif
