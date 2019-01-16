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
class MeshAdapter {
private:
    Mesh& _mesh;

public:
    MeshAdapter(Mesh& mesh) : _mesh(mesh) {}

};

bool try_flip_edge() {
    if(target-already-has-an-edge) {
        issue-a-warning();
        return false;
    }
    if(topo_constraint_not_satisfied) return false;
    if(not-coplanar) return false;
    if(triangle-quality-can-not-be-improved) return false;
    flip_edge();
}

bool try_split_edge() {
    edge_split(); // Bezier curve?
    for_quad_edges_around_the_vertex { try_edge_flip(); }
    return true;
}
bool try_collapse_edge() {
    delta_shape_quality = { try_remove_vertex_1(), try_remove_vertex_2() };
    if(max{delta_shape_quality} < threshold) return false;
    if(deg[0] < deg[1]) {
        remove_vertex_1();
    } else {
        remove_vertex_2();
    }
    return true;
}

void global_relaxation_with_edge_flipping() {
    until( <no-flipping-available> && <relaxation-complete> ) {
        loop_N_times(global_relaxation);
        for_all_edges(try_flipping);
    }
}
void algo() {
    init();

    compute_normal_and_size_measures();
    while( <size-measure-not-satisfied> ) {
        bool topoModified = false;
        do {
            loop_all_edges {
                mark_edge_as_inspection_ready;
            }
            loop_all_edges {
                if( <edge-is-inspection-ready> )
                    if( <edge-too-long> ) {
                        try_split_edge();
                        topoModified = true;
                    } else if( <edge-too-short> ) {
                        try_collapse_edge();
                        topoModified = true;
                    }
            }
        } while( !topoModified );

        global_relaxation_with_edge_flipping();

        compute_normal_and_size_measures();
    };
}
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

    std::transform(
        mesh.getVertices().begin(), mesh.getVertices().end(),
        coords.begin(),
        forces.begin(),
        [&mesh, &coords](const auto& v, const auto& c) {
            mathfunc::Vec3 f {};
            mesh.forEachHalfEdgeTargetingVertex(v, [&](size_t hei) {
                const auto l0 = mesh.getEdgeAttribute(mesh.edge(hei)).adapt.l0;
                const auto r = coords[mesh.target(mesh.opposite(hei))] - c;
                const auto mag = mathfunc::magnitude(r);
                f += r * (1.0 - l0 / mag);
            });

            const auto& un = v.attr.adapt.unitNormal;
            f -= un * dot(un, f); // remove normal component

            return f;
        }
    );

    return forces;
}
struct LocalForceCalculator {
    const std::vector<size_t>& vertexIndex;

    template< typename Mesh >
    auto operator()(const Mesh& mesh, const mathfunc::VecMut< mathfunc::Vec3 >& coords) {
        // The coords must have the same dimension as vertexIndex

        mathfunc::VecMut< mathfunc::Vec3 > forces(coords.size());

        std::transform(
            vertexIndex.begin(), vertexIndex.end(),
            coords.begin(),
            forces.begin(),
            [this, &mesh, &coords](size_t v, const auto& c) {
                mathfunc::Vec3 f {};
                mesh.forEachHalfEdgeTargetingVertex(v, [&](size_t hei) {
                    const auto l0 = mesh.getEdgeAttribute(mesh.edge(hei)).adapt.l0;
                    const auto r = mathfunc::vector2Vec<3, double>(mesh.getVertexAttribute(mesh.target(mesh.opposite(hei))).vertex->coordinate) - c;
                    const auto mag = mathfunc::magnitude(r);
                    f += r * (1.0 - l0 / mag);
                });

                const auto& un = mesh.getVertexAttribute(v).adapt.unitNormal;
                f -= un * dot(un, f); // remove normal component

                return f;
            }
        );

        return forces;
    }
};

auto maxMag2(const mathfunc::VecMut< mathfunc::Vec3 >& v) {
    double res = 0.0;
    for(const auto& i : v) {
        double mag2 = mathfunc::magnitude2(i);
        if(mag2 > res) res = mag2;
    }
    return res;
}

// 2nd order Runge Kutta method on a set of coordinates.
template< typename VecType, typename Func >
void rk2(VecType& coords, Func&& calc_force, double dt, double epsilonSqr) {
    auto forces = calc_force(coords);
    auto maxMag2F = maxMag2(forces);

    while(maxMag2F >= epsilonSqr) {
        // Test move halfway
        auto coords_halfway = coords;
        coords_halfway += forces * (0.5 * dt);

        // Force at halfway
        const auto forces_halfway = calc_force(coords_halfway);

        // Real move
        coords += forces_halfway * dt;

        // Compute new forces
        forces = calc_force(coords);
        maxMag2F = maxMag2(forces);
    }

}

template< typename Mesh > void global_relaxation(Mesh& mesh) {
    using coordinate_type = typename Mesh::VertexAttribute::coordinate_type;

    // TODO Need dt and epsilonSqr
    // normal is obtained from curvature computation
    const size_t numVertices = mesh.numVertices();

    mathfunc::VecMut< mathfunc::Vec3 > coords(numVertices);

    // Copy coordinates to a vector
    std::transform(
        mesh.getVertices().begin(), mesh.getVertices().end(),
        coords.begin(),
        [](const auto& v) { return mathfunc::vector2Vec<3, double>(v.attr.vertex->coordinate); }
    );

    rk2(coords, compute_all_forces, dt, epsilonSqr);

}
template< typename Mesh > void local_relaxation(Mesh& mesh, const std::vector<size_t>& vertexIndex) {
    using coordinate_type = typename Mesh::VertexAttribute::coordinate_type;

    // TODO Need dt and epsilonSqr
    // normal is obtained from curvature computation
    const size_t numVertices = vertexIndex.size();

    mathfunc::VecMut< mathfunc::Vec3 > coords(numVertices);

    // Copy coordinates to a vector
    std::transform(
        vertexIndex.begin(), vertexIndex.end(),
        coords.begin(),
        [&mesh](size_t v) { return mathfunc::vector2Vec<3, double>(mesh.getVertexAttribute(v).vertex->coordinate); }
    );

    rk2(coords, LocalForceCalculator{vertexIndex}, dt, epsilonSqr);

}

// general algorithm procedure
template< typename Mesh >
void adaptive_mesh(Mesh& mesh) {
    calc_l0_all();
    global_relaxation(); // Lowers the force
    edge_flip_list = find_potential_edge_flip();
    for(each_edge_flip : edge_flip_list) do(each_edge_flip);
    loop_through_all_elements(triangles, vertices) {

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
