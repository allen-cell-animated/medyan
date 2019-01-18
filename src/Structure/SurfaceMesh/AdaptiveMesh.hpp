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
#include <cstdint>
#include <limits>
#include <vector>

#include "MathFunctions.h"

#include "Structure/SurfaceMesh/AdaptiveMeshAttribute.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"

enum class TriangleQualityCriteria {
    RadiusRatio     // Circumradius / (2 * Inradius), [1, inf)
};
template< TriangleQualityCriteria > struct TriangleQuality;
template<> struct TriangleQuality< TriangleQualityCriteria::RadiusRatio > {
    static constexpr double best = 1.0;
    static constexpr double worst = std::numeric_limits<double>::infinity();
    static constexpr bool better(double q1, double q2) { return q1 < q2; }
    static constexpr auto betterOne(double q1, double q2) { return better(q1, q2) ? q1 : q2; }
    static constexpr bool worse(double q1, double q2) { return q1 > q2; }
    static constexpr auto worseOne(double q1, double q2) { return worse(q1, q2) ? q1 : q2; }
    static constexpr auto improvement(double q0, double q1) { return q0 / q1; }

    template< typename VecType >
    auto operator()(const VecType& v0, const VecType& v1, const VecType& v2) const {
        using namespace mathfunc;
        const auto d0 = distance(v1, v2);
        const auto d1 = distance(v2, v0);
        const auto d2 = distance(v0, v1);
        const auto p = 0.5 * (d0 + d1 + d2);
        return d0 * d1 * d2 / (8 * (p - d0) * (p - d1) * (p - d2));
    }
}

template< typename Mesh, TriangleQualityCriteria c > class EdgeFlipManager {
public:
    using TriangleQualityType = TriangleQuality< c >;

private:
    size_t _minDegree;
    size_t _maxDegree;
    double _minDotNormal; // To assess whether triangles are coplanar.

public:
    // Returns whether the edge is flipped.
    // Requires
    //   - vertex degrees
    //   - triangle unit normal and quality
    bool tryFlip(Mesh& mesh, size_t ei) const {
        using namespace mathfunc;

        const size_t hei = mesh.getEdges()[ei].halfEdgeIndex;
        const size_t hei_o = mesh.opposite(hei);
        const size_t hei_n = mesh.next(hei);
        const size_t hei_p = mesh.prev(hei);
        const size_t hei_on = mesh.next(hei_o);
        const size_t hei_op = mesh.prev(hei_o);

        const size_t vi0 = mesh.target(hei);
        const size_t vi1 = mesh.target(hei_n);
        const size_t vi2 = mesh.target(hei_o);
        const size_t vi3 = mesh.target(hei_on);
        // Currently the edge connects v0 and v2.
        // If the edge flips, the connection would be between v1 and v3.

        const size_t ti0 = mesh.triangle(hei);
        const size_t ti1 = mesh.triangle(hei_o);

        // Check if topo constraint is satisfied.
        if(
            mesh.degree(vi0) <= _minDegree ||
            mesh.degree(vi2) <= _minDegree ||
            mesh.degree(vi1) >= _maxDegree ||
            mesh.degree(vi3) >= _maxDegree
        ) return false;

        // Check if the current triangles are coplanar.
        if(dot(
            mesh.getTriangleAttribute(ti0).gTriangle.unitNormal,
            mesh.getTriangleAttribute(ti1).gTriangle.unitNormal
        ) < _minDotNormal) return false;

        // Check if the target triangles are coplanar.
        const auto c0 = vector2Vec<3, double>(mesh.getVertexAttribute(vi0).getCoordinate());
        const auto c1 = vector2Vec<3, double>(mesh.getVertexAttribute(vi1).getCoordinate());
        const auto c2 = vector2Vec<3, double>(mesh.getVertexAttribute(vi2).getCoordinate());
        const auto c3 = vector2Vec<3, double>(mesh.getVertexAttribute(vi3).getCoordinate());
        const auto un013 = normalizedVector(cross(c1 - c0, c3 - c0));
        const auto un231 = normalizedVector(cross(c3 - c2, c1 - c2));
        if(dot(un013, un231) < _minDotNormal) return false;

        // Check whether triangle quality can be improved.
        const auto qBefore = TriangleQualityType::betterOne(
            mesh.getTriangleAttribute(ti0).aTriangle.quality,
            mesh.getTriangleAttribute(ti1).aTriangle.quality
        );
        const auto q013 = TriangleQualityType{}(c0, c1, c3);
        const auto q231 = TriangleQualityType{}(c2, c3, c1);
        const auto qAfter = TriangleQualityType::betterOne(q013, q231);
        if( !TriangleQualityType::better(qAfter, qBefore) ) return false;

        // All checks complete. Do the flip.
        typename Mesh::EdgeFlip{}(mesh, ei);

        // TODO: set new triangle attributes
        // TODO: set new edge length (?)
        return true;
    }
};

enum class EdgeSplitVertexInsertionMethod {
    MidPoint
};
template< EdgeSplitVertexInsertionMethod > struct EdgeSplitVertexInsertion;
template<> struct EdgeSplitVertexInsertion< EdgeSplitVertexInsertionMethod::MidPoint > {
    size_t v0, v1;
    template< typename Mesh >
    auto coordinate(const Mesh& mesh, size_t v) const {
        const auto& c0 = mesh.getVertexAttribute(v0).getCoordinate();
        const auto& c1 = mesh.getVertexAttribute(v1).getCoordinate();
        return mathfunc::midPointCoordinate(c0, c1, 0.5);
    }
};

template<
    typename Mesh,
    TriangleQualityCriteria c,
    EdgeSplitVertexInsertionMethod m
> class EdgeSplitManager {
public:
    using EdgeFlipManagerType = EdgeFlipManager< Mesh, c >;
    using EdgeSplitVertexInsertionType = EdgeSplitVertexInsertion< m >;

private:
    size_t _maxDegree;

public:
    // Returns whether a new vertex is inserted.
    // Requires
    //   - Vertex degree
    bool trySplit(Mesh& mesh, size_t ei, const EdgeFlipManagerType& efm) const {
        using namespace mathfunc;

        const size_t hei = mesh.getEdges()[ei].halfEdgeIndex;
        const size_t hei_o = mesh.opposite(hei);
        const size_t hei_n = mesh.next(hei);
        const size_t hei_p = mesh.prev(hei);
        const size_t hei_on = mesh.next(hei_o);
        const size_t hei_op = mesh.prev(hei_o);

        const size_t vi0 = mesh.target(hei);
        const size_t vi1 = mesh.target(hei_n);
        const size_t vi2 = mesh.target(hei_o);
        const size_t vi3 = mesh.target(hei_on);

        const size_t ei0 = mesh.edge(hei_n); // v0 - v1
        const size_t ei1 = mesh.edge(hei_p); // v1 - v2
        const size_t ei2 = mesh.edge(hei_on); // v2 - v3
        const size_t ei3 = mesh.edge(hei_op); // v3 - v1

        // Check topology constraints
        // A new vertex with degree 4 will always be introduced
        if(
            mesh.degree(vi1) >= _maxDegree ||
            mesh.degree(vi3) >= _maxDegree
        ) return false;

        // All checks passed. Do the splitting.
        typename Mesh::VertexInsertionOnEdge< EdgeSplitVertexInsertionType > {}(mesh, ei);

        // TODO: set new triangle attributes
        // TODO: update edge preferred lengths
        // TODO: update edge lengths (?)

        // Propose edge flipping on surrounding quad edges
        efm.tryFlip(mesh, ei0);
        efm.tryFilp(mesh, ei1);
        efm.tryFlip(mesh, ei2);
        efm.tryFlip(mesh, ei3);

        return true;

    }
};

template< typename Mesh, TriangleQualityCriteria c > class EdgeCollapseManager {
public:
    using TriangleQualityType = TriangleQuality< c >;

private:
    size_t _minDegree;
    size_t _maxDegree;
    double _minQualityImprovement; // If smaller than 1, then some degradation is allowed.

public:
    // Returns whether the edge is collapsed
    // Requires
    //   - Triangle quality
    bool tryCollapse(Mesh& mesh, size_t ei) const {
        using namespace mathfunc;

        const size_t hei = mesh.getEdges()[ei].halfEdgeIndex;
        const size_t hei_o = mesh.opposite(hei);
        const size_t hei_n = mesh.next(hei);
        const size_t hei_p = mesh.prev(hei);
        const size_t hei_on = mesh.next(hei_o);
        const size_t hei_op = mesh.prev(hei_o);

        const size_t vi0 = mesh.target(hei);
        const size_t vi1 = mesh.target(hei_n);
        const size_t vi2 = mesh.target(hei_o);
        const size_t vi3 = mesh.target(hei_on);
        // Currently the edge connects v0 and v2.
        // If the edge collapses, v0 and v2 would become one point.

        const size_t ti0 = mesh.triangle(hei);
        const size_t ti1 = mesh.triangle(hei_o);

        // Check topology constraints
        if(
            mesh.degree(vi0) + mesh.degree(vi2) - 4 > _maxDegree ||
            mesh.degree(vi0) + mesh.degree(vi2) - 4 < _minDegree ||
            mesh.degree(vi1) <= _minDegree ||
            mesh.degree(vi3) <= _minDegree
        ) return false;

        // Future: maybe also geometric constraints (gap, smoothness, etc)

        // Check triangle quality constraints
        const auto c0 = vector2Vec<3, double>(mesh.getVertexAttribute(vi0).getCoordinate());
        const auto c2 = vector2Vec<3, double>(mesh.getVertexAttribute(vi2).getCoordinate());

        // Calculate previous triangle qualities around a vertex
        // if v0 is removed
        double q0Before = TriangleQualityType::worst;
        double q0After = TriangleQualityType::worst;
        mesh.forEachHalfEdgeTargetingVertex(vi0, [&](size_t hei) {
            const size_t ti = mesh.triangle(hei);
            q0Before = TriangleQualityType::betterOne(
                mesh.getTriangleAttribute(ti).aTriangle.quality,
                q0Before
            );
            if(ti != ti0 && ti != ti1) {
                const size_t vn = mesh.target(mesh.next(hei));
                const size_t vp = mesh.target(mesh.prev(hei));
                const auto cn = vector2Vec<3, double>(mesh.getVertexAttribute(vn).getCoordinate());
                const auto cp = vector2Vec<3, double>(mesh.getVertexAttribute(vp).getCoordinate());
                q0After = TriangleQualityType::betterOne(
                    TriangleQualityType{}(cp, c2, cn),
                    q0After
                );
            }
        });
        const auto imp0 = TriangleQualityType::improvement(q0Before, q0After);

        // if v2 is removed
        double q2Before = TriangleQualityType::worst;
        double q2After = TriangleQualityType::worst;
        mesh.forEachHalfEdgeTargetingVertex(vi2, [&](size_t hei) {
            const size_t ti = mesh.triangle(hei);
            q2Before = TriangleQualityType::betterOne(
                mesh.getTriangleAttribute(ti).aTriangle.quality,
                q2Before
            );
            if(ti != ti0 && ti != ti1) {
                const size_t vn = mesh.target(mesh.next(hei));
                const size_t vp = mesh.target(mesh.prev(hei));
                const auto cn = vector2Vec<3, double>(mesh.getVertexAttribute(vn).getCoordinate());
                const auto cp = vector2Vec<3, double>(mesh.getVertexAttribute(vp).getCoordinate());
                q2After = TriangleQualityType::betterOne(
                    TriangleQualityType{}(cp, c0, cn),
                    q2After
                );
            }
        });
        const auto imp2 = TriangleQualityType::improvement(q2Before, q2After);

        if(imp0 < _minQualityImprovement && imp2 < _minQualityImprovement) return false;

        if(imp0 > imp2) {
            // Remove v0, collapse onto v2
            typename Mesh::EdgeCollapse {}(mesh, hei_o);
        } else {
            // Remove v2, collapse onto v0
            typename Mesh::EdgeCollapse {}(mesh, hei);
        }

        // TODO: update triangle normals and qualities
        // Do not update edge preferred lengths
        // TODO: update edge lengths?

        return true;
    }
};

template< typename Mesh >
class MeshAdapter {
private:
    Mesh& _mesh;

public:
    MeshAdapter(Mesh& mesh) : _mesh(mesh) {}

};


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
