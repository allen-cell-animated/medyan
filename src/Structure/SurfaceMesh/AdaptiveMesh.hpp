/*

Adaptive mesh algorithm

Implementation inspired by
"An Adaptive Mesh Algorithm for Evolving Surfaces: Simulations of Drop Breakup and Coalescence" (2001)
by Vittorio Cristini, Jerzy Blawzdziewicz and Michael Loewenberg,
"About surface remeshing" (2000) by Pascal J. Frey,
"Geometric surface mesh optimization" (1998) by Pascal J. Frey and Houman Borouchaki

Performs mesh relaxation and topological transformation to improve mesh size
quality and shape quality, while maintaining the geometrical accuracy.

The algorithm works as follows.

Loop
    Update size measures
    Resample vertex to fulfill size quality by inserting/deleting vertex
    Relocate vertex to optimize triangle quality
Until all criteria are met or maximum iterations reached

*/

#ifndef MEDYAN_AdaptiveMesh_hpp
#define MEDYAN_AdaptiveMesh_hpp

#include <algorithm> // max min
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "MathFunctions.h"
#include "Structure/SurfaceMesh/AdaptiveMeshGeometryManager.hpp"
#include "Structure/SurfaceMesh/AdaptiveMeshVertexRelocation.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneMeshCheck.hpp"
#include "Structure/SurfaceMesh/MeshTriangleQuality.hpp"

namespace adaptive_mesh {

// Recommended adaptive mesh parameters
//-------------------------------------
// Topological operations
constexpr size_t surface_mesh_min_degree = 4;
constexpr size_t surface_mesh_max_degree = 9;
constexpr double edge_flip_min_dot_normal = 0.9;
constexpr double edge_collapse_min_quality_improvement = 0.6;
constexpr double edge_collapse_min_dot_normal = 0.85;
// Vertex relocation operations
constexpr double vertex_relaxation_epsilon = 0.05; // (unitless speed/force). The tolerance (l / l_0 - 1)
constexpr double vertex_relaxation_dt = 2.0; // (has unit of length) (around minSize / (iterRelocation * avgForce))
constexpr size_t vertex_relocation_max_iter_relocation = 10;
constexpr size_t vertex_relocation_max_iter_tot = 3; // (vertex relocation + edge flipping) as 1 iter
// Size diffusion
constexpr double size_measure_curvature_resolution = 0.30; // cos of which should be slightly bigger than flip minDotNormal
constexpr double size_measure_max = 60; // Related to the resolution of the system
constexpr size_t size_measure_diffuse_iter = 4;
// Main loop
constexpr size_t mesh_adaptation_topology_max_iter = 8; // Max times of scanning all the edges for sampling adjustment
constexpr size_t mesh_adaptation_soft_max_iter = 8;
constexpr size_t mesh_adaptation_hard_max_iter = 25;

// Implementation
//-------------------------------------

template< typename Mesh, TriangleQualityCriteria c > class EdgeFlipManager {
public:
    using TriangleQualityType = TriangleQuality< c >;
    using CoordinateType = typename Mesh::AttributeType::coordinate_type;

    enum class State {
        Success,
        InvalidTopo,
        NonCoplanar,
        BadQuality
    };

private:
    size_t _minDegree;
    size_t _maxDegree;
    double _minDotNormal; // To assess whether triangles are coplanar.

public:

    // Constructor
    // Parameters
    //   - minDegree: minimum number of neighbors of any vertex
    //   - maxDegree: maximum number of neighbors of any vertex
    //   - minDotNormal: minimum dot product required between neighboring triangles. Range (-1, 1)
    EdgeFlipManager(size_t minDegree, size_t maxDegree, double minDotNormal) :
        _minDegree(minDegree), _maxDegree(maxDegree), _minDotNormal(minDotNormal) {}

    // Returns whether the edge is flipped.
    // Requires
    //   - vertex degrees
    //   - triangle unit normal
    State tryFlip(Mesh& mesh, size_t ei) const {
        using namespace mathfunc;

        const size_t hei = mesh.getEdges()[ei].halfEdgeIndex;
        const size_t hei_o = mesh.opposite(hei);
        const size_t hei_n = mesh.next(hei);
        const size_t hei_on = mesh.next(hei_o);

        const size_t vi0 = mesh.target(hei);
        const size_t vi1 = mesh.target(hei_n);
        const size_t vi2 = mesh.target(hei_o);
        const size_t vi3 = mesh.target(hei_on);
        // Currently the edge connects v0 and v2.
        // If the edge flips, the connection would be between v1 and v3.

        const size_t ti0 = mesh.triangle(hei);
        const size_t ti1 = mesh.triangle(hei_o);

        // Check if topo constraint is satisfied.
        if(mesh.isEdgeOnBorder(ei)) return State::InvalidTopo;

        if(
            mesh.degree(vi0) <= _minDegree ||
            mesh.degree(vi2) <= _minDegree ||
            mesh.degree(vi1) >= _maxDegree ||
            mesh.degree(vi3) >= _maxDegree
        ) return State::InvalidTopo;

        // Check if the current triangles are coplanar.
        if(dot(
            mesh.getTriangleAttribute(ti0).gTriangle.unitNormal,
            mesh.getTriangleAttribute(ti1).gTriangle.unitNormal
        ) < _minDotNormal) return State::NonCoplanar;

        // Check if the target triangles are coplanar.
        const CoordinateType c0 (mesh.getVertexAttribute(vi0).getCoordinate());
        const CoordinateType c1 (mesh.getVertexAttribute(vi1).getCoordinate());
        const CoordinateType c2 (mesh.getVertexAttribute(vi2).getCoordinate());
        const CoordinateType c3 (mesh.getVertexAttribute(vi3).getCoordinate());
        const auto n013 = cross(c1 - c0, c3 - c0);
        const auto mag_n013 = magnitude(n013);
        const auto n231 = cross(c3 - c2, c1 - c2);
        const auto mag_n231 = magnitude(n231);
        if(mag_n013 == 0.0 || mag_n231 == 0.0) return State::BadQuality; // Degenerate case
        if(dot(n013, n231) < _minDotNormal * mag_n013 * mag_n231) return State::NonCoplanar;

        // Check whether triangle quality can be improved.
        const auto q012 = TriangleQualityType{}(c0, c1, c2);
        const auto q230 = TriangleQualityType{}(c2, c3, c0);
        const auto qBefore = TriangleQualityType::worseOne(q012, q230);
        const auto q013 = TriangleQualityType{}(c0, c1, c3);
        const auto q231 = TriangleQualityType{}(c2, c3, c1);
        const auto qAfter = TriangleQualityType::worseOne(q013, q231);
        if( !TriangleQualityType::better(qAfter, qBefore) ) return State::BadQuality;

        // All checks complete. Do the flip.
        typename Mesh::EdgeFlip{}(mesh, ei, [](
            Mesh& mesh, std::array<size_t, 2> tis, std::array<size_t, 4> vis
        ) {
            // Invalidate mesh index cache
            mesh.getMetaAttribute().cacheValid = false;

            for(auto ti : tis) {
                Mesh::AttributeType::adaptiveComputeTriangleNormal(mesh, ti);
                mesh.forEachHalfEdgeInTriangle(ti, [&](size_t hei) {
                    Mesh::AttributeType::adaptiveComputeAngle(mesh, hei);
                });
            }
            for(auto vi : vis) if(!mesh.isVertexOnBorder(vi)) Mesh::AttributeType::adaptiveComputeVertexNormal(mesh, vi);
        });

        // Does not change the edge preferrable length

        return State::Success;
    }
};

enum class EdgeSplitVertexInsertionMethod {
    MidPoint,
    AvgCurv     // Mid point snapped to the sphere of curvature, averaged
};
template< EdgeSplitVertexInsertionMethod > struct EdgeSplitVertexInsertion;
template<> struct EdgeSplitVertexInsertion< EdgeSplitVertexInsertionMethod::MidPoint > {
    size_t v0, v1;
    template< typename Mesh >
    auto coordinate(const Mesh& mesh, size_t v) const {
        using CoordinateType = typename Mesh::AttributeType::coordinate_type;

        const auto c0 = mesh.getVertexAttribute(v0).getCoordinate();
        const auto c1 = mesh.getVertexAttribute(v1).getCoordinate();
        return static_cast<CoordinateType>((c0 + c1) * 0.5);
    }
};
template<> struct EdgeSplitVertexInsertion< EdgeSplitVertexInsertionMethod::AvgCurv > {
    size_t v0, v1;
    // Requires
    //   - Vertex unit normal
    template< typename Mesh >
    auto coordinate(const Mesh& mesh, size_t v) const {
        using namespace mathfunc;
        using CoordinateType = typename Mesh::AttributeType::coordinate_type;

        const CoordinateType c0 (mesh.getVertexAttribute(v0).getCoordinate());
        const CoordinateType c1 (mesh.getVertexAttribute(v1).getCoordinate());
        const auto& un0 = mesh.getVertexAttribute(v0).aVertex.unitNormal;
        const auto& un1 = mesh.getVertexAttribute(v1).aVertex.unitNormal;

        const auto r = c1 - c0;
        const auto mag2_r = magnitude2(r);

        // Compute radius of curvature (for both vertices)
        // negative: convex; positive: concave
        const auto r0 = mag2_r / (2 * dot(un0, r));
        const auto r1 = -mag2_r / (2 * dot(un1, r));

        CoordinateType res0, res1;

        if(std::abs(r0) == std::numeric_limits<double>::infinity()) {
            res0 = 0.5 * (c0 + c1);
        } else {
            // Compute vector from center of sphere to mid point
            const auto ro0 = 0.5 * r - r0 * un0;
            const auto mag_ro0 = magnitude(ro0);

            if(mag_ro0 == 0.0) {
                throw std::runtime_error("Unit normal is parallel to edge.");
            }

            res0 = c0 + r0 * un0 + ro0 * (std::abs(r0) / mag_ro0);
        }

        if(std::abs(r1) == std::numeric_limits<double>::infinity()) {
            res1 = 0.5 * (c0 + c1);
        } else {
            // Compute vector from center of sphere to mid point
            const auto ro1 = -0.5 * r - r1 * un1;
            const auto mag_ro1 = magnitude(ro1);

            if(mag_ro1 == 0.0) {
                throw std::runtime_error("Unit normal is parallel to edge.");
            }

            res1 = c1 + r1 * un1 + ro1 * (std::abs(r1) / mag_ro1);
        }

        return static_cast<CoordinateType>(0.5 * (res0 + res1));
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
    using CoordinateType = typename Mesh::AttributeType::coordinate_type;
    enum class State {
        Success,
        InvalidTopo,
        NotLongestEdge
    };

private:
    size_t _maxDegree;

public:

    // Constructor
    EdgeSplitManager(size_t maxDegree) : _maxDegree(maxDegree) {}

    // Returns whether a new vertex is inserted.
    // Requires
    //   - Vertex degree
    State trySplit(Mesh& mesh, size_t ei, const EdgeFlipManagerType& efm) const {
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
        // Currently does not support insertion on border edges, but we may also implement that in the future.
        if(mesh.isEdgeOnBorder(ei)) return State::InvalidTopo;

        // A new vertex with degree 4 will always be introduced
        if(
            mesh.degree(vi1) >= _maxDegree ||
            mesh.degree(vi3) >= _maxDegree
        ) return State::InvalidTopo;

        // Check whether the current edge is the longest in the triangle
        const CoordinateType c0 (mesh.getVertexAttribute(vi0).getCoordinate());
        const CoordinateType c1 (mesh.getVertexAttribute(vi1).getCoordinate());
        const CoordinateType c2 (mesh.getVertexAttribute(vi2).getCoordinate());
        const CoordinateType c3 (mesh.getVertexAttribute(vi3).getCoordinate());
        const auto l2_e = distance2(c0, c2);
        const auto l2_01 = distance2(c0, c1);
        const auto l2_12 = distance2(c1, c2);
        const auto l2_23 = distance2(c2, c3);
        const auto l2_30 = distance2(c3, c0);
        if( !(
            l2_e >= l2_01 && l2_e >= l2_12 &&
            l2_e >= l2_23 && l2_e >= l2_30
        )) return State::NotLongestEdge;

        // All checks passed. Do the splitting.
        typename Mesh::template VertexInsertionOnEdge< EdgeSplitVertexInsertionType > {}(mesh, ei, [](
            Mesh& mesh, std::array<size_t, 4> tis, std::array<size_t, 5> vis, std::array<size_t, 4> eis
        ) {
            // Invalidate mesh index cache
            mesh.getMetaAttribute().cacheValid = false;

            for(auto ti : tis) {
                Mesh::AttributeType::adaptiveComputeTriangleNormal(mesh, ti);
                mesh.forEachHalfEdgeInTriangle(ti, [&](size_t hei) {
                    Mesh::AttributeType::adaptiveComputeAngle(mesh, hei);
                });
            }
            for(auto vi : vis) if(!mesh.isVertexOnBorder(vi)) Mesh::AttributeType::adaptiveComputeVertexNormal(mesh, vi);

            // Set preferrable length of edges to be the same as before
            const auto eqLength = mesh.getEdgeAttribute(eis[0]).aEdge.eqLength;
            mesh.getEdgeAttribute(eis[1]).aEdge.eqLength = eqLength;
            mesh.getEdgeAttribute(eis[2]).aEdge.eqLength = eqLength;
            mesh.getEdgeAttribute(eis[3]).aEdge.eqLength = eqLength;
        });

        // Propose edge flipping on surrounding quad edges
        efm.tryFlip(mesh, ei0);
        efm.tryFlip(mesh, ei1);
        efm.tryFlip(mesh, ei2);
        efm.tryFlip(mesh, ei3);

        return State::Success;

    }
};

template< typename Mesh, TriangleQualityCriteria c > class EdgeCollapseManager {
public:
    using TriangleQualityType = TriangleQuality< c >;
    using CoordinateType      = typename Mesh::AttributeType::coordinate_type;

    enum class State {
        Success,
        InvalidTopo,
        NonCoplanar,
        BadQuality
    };

private:
    size_t minDegree_;
    size_t maxDegree_;
    double _minQualityImprovement; // If smaller than 1, then some degradation is allowed.
    double _minDotNormal; // Coplanarness requirement after collapse

    struct PrequalifyResult_ {
        double qBefore, qAfter;
        double minCosDihedral;
    };
    // Prequalify the collapse
    // hei is the direction of collapsing (source gets removed, and target is preserved)
    static auto prequalify_(const Mesh& mesh, size_t hei) {
        using namespace mathfunc;

        const auto hei_o = mesh.opposite(hei); // targeting vi1
        const auto hei_n = mesh.next(hei);
        const auto hei_ono = mesh.opposite(mesh.next(hei_o)); // targeting vi1
        const auto hei_opo = mesh.opposite(mesh.prev(hei_o));

        const auto vi0 = mesh.target(hei); // preserved
        const auto vi1 = mesh.target(hei_o); // to be removed
        const CoordinateType c0 (mesh.getVertexAttribute(vi0).getCoordinate());
        const CoordinateType c1 (mesh.getVertexAttribute(vi1).getCoordinate());

        const auto ti0 = mesh.triangle(hei);
        const auto ti1 = mesh.triangle(hei_o);

        double qBefore = TriangleQualityType::best;
        double qAfter = TriangleQualityType::best;
        double minCosDihedral = 1.0; // Coplanar case (best)

        {
            auto chei = hei_o;
            Vec3 lastUnitNormal = mesh.getTriangleAttribute(mesh.triangle(mesh.opposite(hei_n))).gTriangle.unitNormal;
            do {
                const auto ti = mesh.triangle(chei);
                const auto chei_po = mesh.opposite(mesh.prev(chei));
                const auto vn = mesh.target(mesh.next(chei));
                const auto vp = mesh.target(mesh.prev(chei));
                const CoordinateType cn (mesh.getVertexAttribute(vn).getCoordinate());
                const CoordinateType cp (mesh.getVertexAttribute(vp).getCoordinate());

                // Triangle quality before
                qBefore = TriangleQualityType::worseOne(
                    TriangleQualityType{}(cp, c1, cn),
                    qBefore
                );

                // Triangle quality after and Dihedral angle after
                if(ti != ti0 && ti != ti1) {
                    // Quality
                    qAfter = TriangleQualityType::worseOne(
                        TriangleQualityType{}(cp, c0, cn),
                        qAfter
                    );

                    // Dihedral angle
                    auto n_0np = cross(cn - c0, cp - c0);
                    const auto mag_n_0np = magnitude(n_0np);
                    if(mag_n_0np == 0.0) {
                        minCosDihedral = -1.0;
                        break;
                    } else {
                        n_0np *= (1.0 / mag_n_0np);
                        // Dihedral angle with last triangle
                        minCosDihedral = std::min(
                            minCosDihedral,
                            dot(n_0np, lastUnitNormal)
                        );
                        lastUnitNormal = n_0np;
                        // Dihedral angle with outside triangle
                        minCosDihedral = std::min(
                            minCosDihedral,
                            dot(n_0np, mesh.getTriangleAttribute(mesh.triangle(chei_po)).gTriangle.unitNormal)
                        );
                        // Special dihedral angle
                        if(chei == hei_ono) {
                            minCosDihedral = std::min(
                                minCosDihedral,
                                dot(n_0np, mesh.getTriangleAttribute(mesh.triangle(hei_opo)).gTriangle.unitNormal)
                            );
                        }
                    }
                }

                // Change chei
                chei = mesh.prev(mesh.opposite(chei)); // counter clockwise around vi1
            } while(chei != hei_o);
        }

        return PrequalifyResult_{ qBefore, qAfter, minCosDihedral };
    }

public:

    // Constructor
    EdgeCollapseManager(size_t minDegree, size_t maxDegree, double minQualityImprovement, double minDotNormal) :
        minDegree_(minDegree), maxDegree_(maxDegree), _minQualityImprovement(minQualityImprovement), _minDotNormal(minDotNormal) {}

    // Returns whether the edge is collapsed
    // Requires
    //   - <None>
    State tryCollapse(Mesh& mesh, size_t ei) const {
        using namespace mathfunc;

        const size_t hei = mesh.getEdges()[ei].halfEdgeIndex;
        const size_t hei_o = mesh.opposite(hei);
        const size_t hei_n = mesh.next(hei);
        const size_t hei_on = mesh.next(hei_o);

        const size_t vi0 = mesh.target(hei);
        const size_t vi1 = mesh.target(hei_n);
        const size_t vi2 = mesh.target(hei_o);
        const size_t vi3 = mesh.target(hei_on);
        // Currently the edge connects v0 and v2.
        // If the edge collapses, v0 and v2 would become one point.

        // Check topology constraints
        // Currently does not allow collapsing of border edges, but we may also implement that in the future
        if(mesh.isEdgeOnBorder(ei)) return State::InvalidTopo;

        if(
            mesh.degree(vi0) + mesh.degree(vi2) - 4 > maxDegree_ ||
            mesh.degree(vi0) + mesh.degree(vi2) - 4 < minDegree_ ||
            (mesh.isVertexOnBorder(vi0) && mesh.isVertexOnBorder(vi2)) ||
            mesh.degree(vi1) <= minDegree_ ||
            mesh.degree(vi3) <= minDegree_
        ) return State::InvalidTopo;

        // Check triangle quality constraints
        // Calculate previous triangle qualities around a vertex
        // if v0 is removed
        const auto prequalifyResult0 = prequalify_(mesh, hei_o);
        const double q0Before = prequalifyResult0.qBefore;
        const double q0After  = prequalifyResult0.qAfter;
        const auto imp0 = TriangleQualityType::improvement(q0Before, q0After);
        const double minCosDihedral0 = prequalifyResult0.minCosDihedral;

        // if v2 is removed
        const auto prequalifyResult2 = prequalify_(mesh, hei);
        const double q2Before = prequalifyResult2.qBefore;
        const double q2After  = prequalifyResult2.qAfter;
        const auto imp2 = TriangleQualityType::improvement(q2Before, q2After);
        const double minCosDihedral2 = prequalifyResult2.minCosDihedral;

        if(imp0 < _minQualityImprovement && imp2 < _minQualityImprovement) return State::BadQuality;
        if(
            (imp0 > imp2 && minCosDihedral0 < _minDotNormal) ||
            (imp0 <= imp2 && minCosDihedral2 < _minDotNormal)
        ) return State::NonCoplanar;

        auto attributeSetter = [](
            Mesh& mesh, size_t hei_begin, size_t hei_end, size_t ov0
        ) {
            // Invalidate mesh index cache
            mesh.getMetaAttribute().cacheValid = false;

            for(size_t hei1 = hei_begin; hei1 != hei_end; hei1 = mesh.opposite(mesh.next(hei1))) {
                const size_t hei1_n = mesh.next(hei1);
                const size_t ti = mesh.triangle(hei1);
                Mesh::AttributeType::adaptiveComputeTriangleNormal(mesh, ti);
                mesh.forEachHalfEdgeInTriangle(ti, [&](size_t hei) {
                    Mesh::AttributeType::adaptiveComputeAngle(mesh, hei);
                });
            }

            if(!mesh.isVertexOnBorder(ov0)) Mesh::AttributeType::adaptiveComputeVertexNormal(mesh, ov0);
            const auto v_first = mesh.target(mesh.opposite(hei_begin));
            if(!mesh.isVertexOnBorder(v_first)) Mesh::AttributeType::adaptiveComputeVertexNormal(mesh, v_first);
            for(size_t hei1 = hei_begin; hei1 != hei_end; hei1 = mesh.opposite(mesh.next(hei1))) {
                const auto vi = mesh.target(mesh.next(hei1));
                if(!mesh.isVertexOnBorder(vi))
                    Mesh::AttributeType::adaptiveComputeVertexNormal(mesh, vi);
            }
        };

        if(imp0 > imp2) {
            // Remove v0, collapse onto v2
            typename Mesh::EdgeCollapse {}(mesh, hei_o, attributeSetter);
        } else {
            // Remove v2, collapse onto v0
            typename Mesh::EdgeCollapse {}(mesh, hei, attributeSetter);
        }

        // Does not update edge preferred lengths

        return State::Success;
    }
};

enum class SizeMeasureCriteria {
    Curvature
};
template< SizeMeasureCriteria > struct VertexSizeMeasure;
template<> struct VertexSizeMeasure< SizeMeasureCriteria::Curvature > {
    double resolution; // size = res * min_radius_curvature
    double upperLimit; // maximum size

    // Requires
    //   - Vertex unit normal
    template< typename Mesh > auto vertexSize(Mesh& mesh, size_t vi) const {
        using CoordinateType = typename Mesh::AttributeType::coordinate_type;

        double minRadiusCurvature = std::numeric_limits<double>::infinity();
        const auto& un = mesh.getVertexAttribute(vi).aVertex.unitNormal;
        const CoordinateType ci (mesh.getVertexAttribute(vi).getCoordinate());
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](size_t hei) {
            const auto r = mesh.getVertexAttribute(mesh.target(mesh.opposite(hei))).getCoordinate() - ci;
            minRadiusCurvature = std::min(
                std::abs(0.5 * mathfunc::magnitude2(r) / mathfunc::dot(un, r)),
                minRadiusCurvature
            );
        });
        
        return std::min(resolution * minRadiusCurvature, upperLimit);
    }
};

template< SizeMeasureCriteria... > struct VertexSizeMeasureCombined;
template< SizeMeasureCriteria c, SizeMeasureCriteria... cs >
struct VertexSizeMeasureCombined< c, cs... > {
    template< typename Mesh >
    static auto vertexSize(Mesh& mesh, size_t vi, const VertexSizeMeasure<c>& vsm, const VertexSizeMeasure<cs>&... vsms) {
        return std::min(vsm.vertexSize(mesh, vi), VertexSizeMeasureCombined<cs...>::vertexSize(mesh, vi, vsms...));
    }
};
template< SizeMeasureCriteria c >
struct VertexSizeMeasureCombined< c > {
    template< typename Mesh >
    static auto vertexSize(Mesh& mesh, size_t vi, const VertexSizeMeasure<c>& vsm) {
        return vsm.vertexSize(mesh, vi);
    }
};

template< typename Mesh > class SizeMeasureManager {
private:

    double _curvRes; // resolution used in radius curvature
    double _maxSize; // Hard upper bound of size
    size_t _diffuseIter; // Diffusion iterations used in gradation control

    template< SizeMeasureCriteria... cs >
    auto _vertexSize(Mesh& mesh, size_t vi, const VertexSizeMeasure<cs>&... vsms) const {
        return VertexSizeMeasureCombined<cs...>::vertexSize(mesh, vi, vsms...);
    }
    template< SizeMeasureCriteria... cs >
    void _updateVertexSize(Mesh& mesh, const VertexSizeMeasure<cs>&... vsms) const {
        const size_t numVertices = mesh.getVertices().size();
        for(size_t i = 0; i < numVertices; ++i) {
            mesh.getVertexAttribute(i).aVertex.size = _vertexSize(mesh, i, vsms...);
        }
    }

    void _diffuseSize(Mesh& mesh) const {
        const size_t numVertices = mesh.getVertices().size();

        // Initialize with max size
        for(size_t i = 0; i < numVertices; ++i) {
            auto& av = mesh.getVertexAttribute(i).aVertex;
            av.size = std::min(av.size, _maxSize);
        }

        // Diffuse, with D * Delta t = 0.5, and uniformly weighted Laplace operator
        // l_new = l_old / 2 + (sum of neighbor l_old) / (2 * numNeighbors)
        for(size_t iter = 0; iter < _diffuseIter; ++iter) {
            for(size_t i = 0; i < numVertices; ++i) {
                auto& av = mesh.getVertexAttribute(i).aVertex;
                const size_t deg = mesh.degree(i);
    
                double sumSizeNeighbor = 0.0;
                mesh.forEachHalfEdgeTargetingVertex(i, [&](size_t hei) {
                    sumSizeNeighbor += mesh.getVertexAttribute(mesh.target(mesh.opposite(hei))).aVertex.size;
                });

                av.sizeAux = std::min(
                    0.5 * av.size + 0.5 * sumSizeNeighbor / deg,
                    _maxSize
                ); // capped by _maxSize
            }
            for(size_t i = 0; i < numVertices; ++i) {
                auto& av = mesh.getVertexAttribute(i).aVertex;
                av.size = av.sizeAux;
            }
        }
    }

    void _updateEdgeEqLength(Mesh& mesh) const {
        const size_t numEdges = mesh.getEdges().size();
        for(size_t i = 0; i < numEdges; ++i) {
            auto& l0 = mesh.getEdgeAttribute(i).aEdge.eqLength;
            l0 = 0.0;
            mesh.forEachHalfEdgeInEdge(i, [&](size_t hei) {
                l0 += 0.5 * mesh.getVertexAttribute(mesh.target(hei)).aVertex.size;
            });
        }
    }

public:

    // Constructor
    SizeMeasureManager(double curvRes, double maxSize, size_t diffuseIter) :
        _curvRes(curvRes), _maxSize(maxSize), _diffuseIter(diffuseIter) {}

    // Requires
    //   - Unit normal on each vertex
    void computeSizeMeasure(Mesh& mesh) const {
        VertexSizeMeasure< SizeMeasureCriteria::Curvature > vsmCurv {_curvRes, _maxSize};

        // Compute size on each vertex
        _updateVertexSize(mesh, vsmCurv);

        // Diffuse size on vertices
        _diffuseSize(mesh);

        // Compute preferred length of edges
        _updateEdgeEqLength(mesh);
    }
}; // End SizeMesasureManager

template< typename Mesh >
class MeshAdapter {
public:
    using GeometryManagerType = GeometryManager< Mesh >;
    using CoordinateType      = typename Mesh::AttributeType::coordinate_type;

    static constexpr auto optimalVertexLocationMethod    = OptimalVertexLocationMethod::barycenter;
    static constexpr auto triangleQualityCriteria        = TriangleQualityCriteria::radiusRatio;
    static constexpr auto edgeSplitVertexInsertionMethod = EdgeSplitVertexInsertionMethod::AvgCurv;

    struct Parameter {
        // Topology
        size_t minDegree;
        size_t maxDegree;
        double edgeFlipMinDotNormal;
        double edgeCollapseMinQualityImprovement;
        double edgeCollapseMinDotNormal;

        // Relaxation
        double relaxationEpsilon;
        double relaxationDt;
        size_t relaxationMaxIterRelocation;
        size_t relaxationMaxIterTotal;

        // Size diffusion
        double curvatureResolution;
        double maxSize;
        size_t diffuseIter;

        // Main loop
        size_t samplingAdjustmentMaxIter;
        size_t mainLoopSoftMaxIter;
        size_t mainLoopHardMaxIter;
    };

private:
    SizeMeasureManager< Mesh > _sizeMeasureManager;
    DirectVertexRelocationManager< Mesh, optimalVertexLocationMethod > _directVertexRelocationManager;

    EdgeFlipManager< Mesh, triangleQualityCriteria > _edgeFlipManager;
    EdgeSplitManager< Mesh, triangleQualityCriteria, edgeSplitVertexInsertionMethod > _edgeSplitManager;
    EdgeCollapseManager< Mesh, triangleQualityCriteria > _edgeCollapseManager;

    size_t _samplingAdjustmentMaxIter; // Maximum number of scans used in sampling.
    size_t _mainLoopSoftMaxIter; // Maximum iterations of the main loop if topology changes can be reduced to 0
    size_t _mainLoopHardMaxIter; // Maximum iterations of the main loop (hard cap)

    void computeSizeMeasures_(Mesh& mesh) const {
        GeometryManagerType::computeAllTriangleNormals(mesh);
        GeometryManagerType::computeAllAngles(mesh);
        GeometryManagerType::computeAllVertexNormals(mesh);
        _sizeMeasureManager.computeSizeMeasure(mesh);
    }
public:
    // Constructor
    MeshAdapter(Parameter param) :
        _sizeMeasureManager(param.curvatureResolution, param.maxSize, param.diffuseIter),
        _directVertexRelocationManager(
            param.relaxationMaxIterRelocation,
            param.relaxationMaxIterTotal
        ),
        _edgeFlipManager(param.minDegree, param.maxDegree, param.edgeFlipMinDotNormal),
        _edgeSplitManager(param.maxDegree),
        _edgeCollapseManager(
            param.minDegree,
            param.maxDegree,
            param.edgeCollapseMinQualityImprovement,
            param.edgeCollapseMinDotNormal
        ),
        _samplingAdjustmentMaxIter(param.samplingAdjustmentMaxIter),
        _mainLoopSoftMaxIter(param.mainLoopSoftMaxIter),
        _mainLoopHardMaxIter(param.mainLoopHardMaxIter)
    {}

    void adapt(Mesh& mesh) const {
        using namespace mathfunc;

        size_t mainLoopIter = 0;
        while(true) {
            // This outer loop does the following in sequence:
            //
            // 1. Update the desired sizes of all the edges, based on the
            //    size criteria.
            // 2. Perform vertex insertion/deletion operations according to the
            //    desired size. The geometry should be updated.
            // 3. Relocate the vertices to local optimum iteratively. The
            //    geometry should be updated.

            computeSizeMeasures_(mesh);

            bool sizeMeasureSatisfied = true;

            size_t countTopoModified;
            size_t iter = 0;
            do {
                // The inner loop traverses the meshwork multiple times to
                // check for elements that do not satisfy the size criteria,
                // and try to fix them by vertex insertion/deletion operations.

                countTopoModified = 0;
                for(size_t ei = 0; ei < mesh.getEdges().size(); /* No increment here */) {
                    const size_t hei0 = mesh.getEdges()[ei].halfEdgeIndex;
                    const size_t v0 = mesh.target(hei0);
                    const size_t v1 = mesh.target(mesh.opposite(hei0));

                    const CoordinateType c0 (mesh.getVertexAttribute(v0).getCoordinate());
                    const CoordinateType c1 (mesh.getVertexAttribute(v1).getCoordinate());
                    const double length2 = distance2(c0, c1);

                    const double eqLength = mesh.getEdgeAttribute(ei).aEdge.eqLength;
                    const double eqLength2 = eqLength * eqLength;

                    if(length2 >= 2 * eqLength2) { // Too long
                        sizeMeasureSatisfied = false;
                        if(_edgeSplitManager.trySplit(mesh, ei, _edgeFlipManager) == decltype(_edgeSplitManager)::State::Success)
                            // Edge splitting happened. Will check edge ei again next round
                            ++countTopoModified;
                        else
                            ++ei;
                    } else if(2 * length2 <= eqLength2) { // Too short
                        sizeMeasureSatisfied = false;
                        if(_edgeCollapseManager.tryCollapse(mesh, ei) == decltype(_edgeCollapseManager)::State::Success)
                            // Edge collapsing happened. The edge at ei will be different next round
                            ++countTopoModified;
                        else
                            ++ei;
                    } else { // Check passed
                        ++ei;
                    }
                }

                ++iter;
            } while(countTopoModified && iter < _samplingAdjustmentMaxIter); // If any topology was modified, will loop through all edges again.

            if(
                sizeMeasureSatisfied
                || (
                    countTopoModified == 0
                    && mainLoopIter >= _mainLoopSoftMaxIter
                )
                || mainLoopIter >= _mainLoopHardMaxIter
            ) break;

            _directVertexRelocationManager(mesh, _edgeFlipManager);

            ++mainLoopIter;
        } // End loop TopoModifying-Relaxation

    } // End function adapt(...)

};

using MembraneMeshAdapter = MeshAdapter< typename Membrane::MeshType >;

} // namespace adaptive_mesh

#endif
