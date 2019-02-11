#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshCheck_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshCheck_hpp

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "util/io/log.h"

namespace membrane_mesh_check {

using MeshType = Membrane::MeshType;

struct MembraneTopologyCheck {
    size_t minDegree;
    size_t maxDegree;
    size_t genus = 0;
    size_t numBoundaries = 0;

    bool operator()(const MeshType& mesh, bool report = false) const {
        bool res = true;

        const size_t chi = 2 - 2 * genus - numBoundaries; // Euler characteristic

        // Check number of elements
        const size_t numVertices = mesh.getVertices().size();
        const size_t numTriangles = mesh.getTriangles().size();
        const size_t numEdges = mesh.getEdges().size();
        const size_t numHalfEdges = mesh.getHalfEdges().size();

        const size_t chiActual = numTriangles + numVertices - numEdges;
        if(chiActual != chi) {
            res = false;
            if(report) {
                LOG(ERROR) << "Incorrect Euler characteristic (expected " << chi << "): " << chiActual;
            }
        }
        if(!numBoundaries && (numEdges * 2 != numTriangles * 3 || numHalfEdges != numEdges * 2)) {
            res = false;
            if(report) {
                LOG(ERROR) << "Incorrect number of elements of a closed surface";
            }
        }

        // Check vertices
        for(size_t i = 0; i < numVertices; ++i) {
            // Check degree
            const size_t degree = mesh.degree(i);
            size_t degreeActual = 0;
            mesh.forEachHalfEdgeTargetingVertex(i, [&](size_t hei) {
                ++degreeActual;
            });
            if(degreeActual != degree) {
                res = false;
                if(report)
                    LOG(ERROR) << "Inconsistent degree at vertex " << i << ": " << degreeActual
                        << " (expected " << degree << ")";
            }
            if(degree < minDegree) {
                res = false;
                if(report)
                    LOG(ERROR) << "Vertex " << i << " has too small degree " << degree;
            }
            if(degree > maxDegree) {
                res = false;
                if(report)
                    LOG(ERROR) << "Vertex " << i << " has too large degree " << degree;
            }
            // Check targeting half edge index
            const size_t hei = mesh.getVertices()[i].halfEdgeIndex;
            const size_t vi = mesh.target(hei);
            if(vi != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent targeting half edge at vertex " << i;
                    LOG(INFO) << "Half edge " << hei << " points to vertex " << vi;
                }
            }
        }

        // Check edges
        for(size_t i = 0; i < numEdges; ++i) {
            // Check half edge index
            const size_t hei = mesh.getEdges()[i].halfEdgeIndex;
            const size_t ei = mesh.edge(hei);
            if(ei != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent half edge index at edge " << i;
                    LOG(INFO) << "Half edge " << hei << " has edge " << ei;
                }
            }
        }

        // Check triangles
        for(size_t i = 0; i < numTriangles; ++i) {
            // Check half edge index
            const size_t hei = mesh.getTriangles()[i].halfEdgeIndex;
            const size_t ti = mesh.triangle(hei);
            if(ti != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent half edge index at triangle " << i;
                    LOG(INFO) << "Half edge " << hei << " has triangle " << ti;
                }
            }
        }

        // Check half edges
        for(size_t i = 0; i < numHalfEdges; ++i) {
            // Check opposite
            if(mesh.hasOpposite(i)) {
                const size_t hei_o = mesh.opposite(i);
                if(!mesh.hasOpposite(hei_o)) {
                    res = false;
                    if(report) {
                        LOG(ERROR) << "Half edge " << hei_o << " (opposite of half edge " << i << ") "
                            << "does not have a half edge.";
                    }
                } else if(mesh.opposite(hei_o) != i) {
                    res = false;
                    if(report) {
                        LOG(ERROR) << "Inconsistent opposite half edges: "
                            << i << " -> " << hei_o << " -> " << mesh.opposite(hei_o);
                    }
                }
            } else {
                // Does not have an opposite
                if(!numBoundaries) {
                    res = false;
                    if(report) {
                        LOG(ERROR) << "In closed surface, half edge " << i << " does not have an opposite.";
                    }
                }
            }
            // Check next/prev
            const size_t next = mesh.next(i);
            const size_t prev = mesh.prev(i);
            if(i != mesh.prev(next)) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Next of half edge " << i << " is " << next
                        << ", but its prev is " << mesh.prev(next);
                }
            }
            if(i != mesh.next(prev))  {
                res = false;
                if(report) {
                    LOG(ERROR) << "Prev of half edge " << i << " is " << prev
                        << ", but its next is " << mesh.next(prev);
                }
            }
        }

        if(!res && report)
            LOG(INFO)
                << "Triangles: " << numTriangles
                << "Vertices: " << numVertices
                << "Edges: " << numEdges
                << "Half edges: " << numHalfEdges;
    }

};

} // namespace membrane_mesh_check

#endif
