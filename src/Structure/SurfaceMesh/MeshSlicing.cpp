#include "MeshSlicing.h"

#include "SimpleGeometricObject.h"
#include "Membrane.h"
#include "Triangle.h"
#include "Vertex.h"

#include "MathFunctions.h"
using namespace mathfunc;

void MeshSlicingManager::planeSliceTriangle(size_t aspect, double otherCoord, Triangle* triangle) {
    /**************************************************************************
    After slicing, what will be produced:
        - Triangles will be sliced and contain 0 or more more 3d polygons
    
    The calculation requires
        - The unit normal of the interested triangle
    
    For simplicity, currently slicing is only valid for CONVEX polygons.
    **************************************************************************/

    // Check whether this triangle can be cut
    auto& v0 = triangle->getVertices()[0]->coordinate;
    auto& v1 = triangle->getVertices()[1]->coordinate;
    auto& v2 = triangle->getVertices()[2]->coordinate;
    if((v0[aspect] - otherCoord) * (v1[aspect] - otherCoord) > 0
        && (v1[aspect] - otherCoord) * (v2[aspect] - otherCoord) > 0)
        return; // No need to cut this triangle
    
    // Get normal vector
    auto& normal = triangle->getGTriangle()->getUnitNormal();
    
    // Check whether this triangle has already been cut
    auto& po = triangle->getGTriangle()->getPolygons();
    if(po.size() == 0) {
        // Create the first polygon using this triangle
        auto sv0 = make_unique<SimpleVertex<3>>(vector2Array<double, 3>(v0));
        auto sv1 = make_unique<SimpleVertex<3>>(vector2Array<double, 3>(v1));
        auto sv2 = make_unique<SimpleVertex<3>>(vector2Array<double, 3>(v2));
        po.emplace_back(vector<SimpleVertex<3>*>{{sv0.get(), sv1.get(), sv2.get()}});
        po.back().setUnitNormal(normal);

        vertices3.insert(move(sv0));
        vertices3.insert(move(sv1));
        vertices3.insert(move(sv2));
    }

    size_t numPolygons = po.size();
    for(size_t idx = 0; idx < numPolygons; ++idx) {
        // Note that now po.size() could change, but the new ones are added at the end of the vector

        double firstSign = po[idx].vertex(0)[aspect] - otherCoord;
        double lastSign = firstSign;
        size_t vIdx = 0;
        SimplePolygon<3>* lastNewPolygon = nullptr;
        // Search and find 0 or 2 intersections. Does not work for concave polygons
        while(vIdx < po[idx].getVertexNum()) {
            double nextSign = po[idx].vertex((vIdx + 1) % po[idx].getVertexNum())[aspect] - otherCoord;

            // Check whether an intersection occurs
            if(lastSign * nextSign < 0) {
                // An intersection occurs

                // Find coordinate of intersection
                double lastShare = lastSign / (lastSign - nextSign);
                double nextShare = nextSign / (nextSign - lastSign);
                array<double, 3> vIntersect;
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    if(coordIdx == aspect) {
                        vIntersect[coordIdx] = otherCoord;
                    } else {
                        vIntersect[coordIdx] = lastShare * po[idx].vertex(vIdx)[coordIdx]
                            + nextShare * po[idx].vertex((vIdx + 1) % po[idx].getVertexNum())[coordIdx];
                    }
                }

                // New simple vertex at intersection
                auto sv = make_unique<SimpleVertex<3>>(vIntersect);

                // Check whether this marks the start or the end of the new polygon
                if(lastNewPolygon) {
                    // This is the end of the new polygon

                    // Give the vertex to new polygon
                    lastNewPolygon->addVertex(&(po[idx].vertex(vIdx)), lastNewPolygon->getVertexNum());
                    po[idx].removeVertex(vIdx);

                    // Insert the new vertex for both polygons
                    lastNewPolygon->addVertex(sv.get(), lastNewPolygon->getVertexNum());
                    po[idx].addVertex(sv.get(), vIdx);

                    // Untrack last polygon
                    lastNewPolygon = nullptr;

                    ++vIdx;

                } else {
                    // This is the start of the new polygon

                    // Add an empty polygon
                    po.emplace_back();
                    lastNewPolygon = &(po.back());
                    lastNewPolygon->setUnitNormal(normal);

                    // Insert the new vertex for both polygons
                    lastNewPolygon->addVertex(sv.get());
                    po[idx].addVertex(sv.get(), vIdx + 1);

                    ++vIdx;
                }

                // Store the new vertex
                vertices3.insert(move(sv));

            } else {
                // No intersection occurs

                if(lastNewPolygon) {
                    // Give the vertex to new polygon
                    lastNewPolygon->addVertex(&(po[idx].vertex(vIdx)), lastNewPolygon->getVertexNum());
                    po[idx].removeVertex(vIdx);

                    // Note that vIdx do not increase in this case!
                } else {
                    ++vIdx;
                }
            }

            lastSign = nextSign;

        } // End of loop on the vertices of the polygon
    } // End of loop on several polygons
}

PlaneSliceSnapshot MeshSlicingManager::planeSliceMembrane(size_t aspect, double otherCoord, unordered_set<Triangle*>& triangles) {
    /**************************************************************************
    This function produces a cut snapshot by iteratively finding triangle cuts
    from the preset set. Which triangles are needed are decided by the caller
    of this function.

    If a cut stopped without forming a complete polygon, an exception will be
    thrown.

    Note that this function is destructive to the original "triangles" set.
    **************************************************************************/
}
