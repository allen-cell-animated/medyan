#include "MeshSlicing.h"

#include <stdexcept>

#include "SimpleGeometricObject.h"
#include "Membrane.h"
#include "GTriangle.h"
#include "Triangle.h"
#include "Vertex.h"
#include "Edge.h"

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
    from the preset set. The triangles needed for search are decided by the
    caller of this function.

    Note that this function is destructive to the passed "triangles" set.

    The function requires:
        - The unit normal of the triangles
    **************************************************************************/

    PlaneSliceSnapshot res(aspect, otherCoord);

    while(!triangles.empty()) {
        unordered_set<Triangle*>::iterator curTriangle = triangles.begin();
        SimplePolygon<2>* np = nullptr;

        Edge* firstEdge = nullptr;
        Edge* nextEdge = nullptr;
        Triangle* lastTriangle = nullptr;
        // Deal with the first triangle
        bool oneTriangleSliced = false;
        while(!oneTriangleSliced) {
            auto& v0 = (*curTriangle)->getVertices()[0]->coordinate;
            auto& v1 = (*curTriangle)->getVertices()[1]->coordinate;
            auto& v2 = (*curTriangle)->getVertices()[2]->coordinate;

            double d0 = v0[aspect] - otherCoord;
            double d1 = v1[aspect] - otherCoord;
            double d2 = v2[aspect] - otherCoord;

            if(d0 * d1 > 0 && d1 * d2 > 0) {
                // This triangle cannot be sliced
                triangles.erase(curTriangle);

                if(triangles.empty()) {
                    // No more triangles to be searched
                    return res;
                }

            } else {
                // This triangle can be sliced
                
                // Create a new polygon
                res.polygons.emplace_back();
                np = &(res.polygons.back());

                // Find the place of the cut
                vector<size_t> edgeId;
                if(d0 * d1 <= 0) edgeId.push_back(0);
                if(d1 * d2 <= 0) edgeId.push_back(1);
                if(d2 * d0 <= 0) edgeId.push_back(2);

                if(edgeId.size() != 2) throw std::runtime_error("The plane must intersect the triangle on exactly 2 edges.");

                for(size_t eachEdgeId: edgeId) {
                    array<double, 2> vIntersect;
                    for(size_t coordIdx = 0; coordIdx < 2; ++coordIdx) {
                        size_t actualCoordIdx = (coordIdx + 1 + aspect) % 3;
                        switch(eachEdgeId) {
                        case 0:
                            vIntersect[coordIdx] = d0 / (d0 - d1) * v0[actualCoordIdx] + d1 / (d1 - d0) * v1[actualCoordIdx];
                            break;
                        case 1:
                            vIntersect[coordIdx] = d1 / (d1 - d2) * v1[actualCoordIdx] + d2 / (d2 - d1) * v2[actualCoordIdx];
                            break;
                        case 2:
                            vIntersect[coordIdx] = d2 / (d2 - d0) * v2[actualCoordIdx] + d0 / (d0 - d2) * v0[actualCoordIdx];
                            break;
                        }
                    }

                    // New simple vertex at intersection
                    auto sv = make_unique<SimpleVertex<2>>(vIntersect);
                    np->addVertex(sv.get(), np->getVertexNum());
                    vertices2.insert(move(sv));
                }

                auto vectorZeroToOne = vectorDifference(np->vertex(1)->getCoordinate(), np->vertex(0)->getCoordinate());
                auto& triangleNormal = (*curTriangle)->getGTriangle()->getUnitNormal();

                if(dotProduct(crossProduct(triangleNormal, vectorZeroToOne), planeNormal[aspect]) < 0) {
                    // Swap vertex 0 and 1 in polygon
                    swap(np->vertex(0), np->vertex(1));
                    firstEdge = (*curTriangle)->getEdges()[edgeId[1]];
                    nextEdge = (*curTriangle)->getEdges()[edgeId[0]];
                } else {
                    firstEdge = (*curTriangle)->getEdges()[edgeId[0]];
                    nextEdge = (*curTriangle)->getEdges()[edgeId[1]];
                }

                lastTriangle = *curTriangle;
                triangles.erase(curTriangle);

                // Mark as already sliced
                oneTriangleSliced = true;
            }           
        }

        // Till this place, we should have non-nullptr value for np, firstEdge, nextEdge, lastTriangle

        bool completePolygon = false;

        while(!completePolygon) {

            // Find the current triangle
            Triangle* nowTriangle = nullptr;
            for(auto eachTriangle: nextEdge->getTriangles()) {
                if(eachTriangle != lastTriangle) {
                    nowTriangle = eachTriangle;
                    break;
                }
            }
            
            // Find which edge is responsible for new intersection
            for(Edge* eachEdge: nowTriangle->getEdges()) {
                if(eachEdge != nextEdge) {
                    auto& v0 = eachEdge->getVertices()[0]->coordinate;
                    auto& v1 = eachEdge->getVertices()[1]->coordinate;
                    double d0 = v0[aspect] - otherCoord;
                    double d1 = v1[aspect] - otherCoord;
                    if(d0 * d1 <= 0) { // Cut on this edge

                        // Check whether the new edge is the first edge
                        if(eachEdge == firstEdge) {
                            completePolygon = true;
                            // Do not add new vertices
                        } else {
                            // Find the intersect and add the vertex to polygon

                            array<double, 2> vIntersect;
                            for(size_t coordIdx = 0; coordIdx < 2; ++coordIdx) {
                                size_t actualCoordIdx = (coordIdx + 1 + aspect) % 3;
                                vIntersect[coordIdx] = d0 / (d0 - d1) * v0[actualCoordIdx] + d1 / (d1 - d0) * v1[actualCoordIdx];
                            }

                            // New simple vertex at intersection
                            auto sv = make_unique<SimpleVertex<2>>(vIntersect);
                            np->addVertex(sv.get(), np->getVertexNum());
                            vertices2.insert(move(sv));

                            // Clear this triangle from the set
                            unordered_set<Triangle*>::iterator itTriangle = triangles.find(nowTriangle);
                            if(itTriangle != triangles.end()) triangles.erase(itTriangle);

                        }

                        // Assign the edge and triangle and exit
                        nextEdge = eachEdge;
                        lastTriangle = nowTriangle;
                        break;
                    }
                }
            }
        } // End of loop of adding vertices on the polygon
    }

    return res;
}
