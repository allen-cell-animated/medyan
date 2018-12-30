#include "Membrane.h"

#include <limits>
#include <stdexcept>
#include <unordered_set>

#include "common.h"
#include "MathFunctions.h"

#include "Compartment.h"
#include "core/controller/GController.h"
#include "SubSystem.h"
#include "SysParams.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"
#include "GTriangle.h"
#include "MTriangle.h"
#include "GEdge.h"
#include "GVoronoiCell.h"
#include "MVoronoiCell.h"
#include "GMembrane.h"
#include "MMembrane.h"

using namespace mathfunc;

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(
    SubSystem* s,
    short membraneType,
    const std::vector< MembraneMeshAttribute::coordinate_type >& vertexCoordinateList,
    const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
) : Trackable(false, false, false, false, true), Geometric(),
    _mesh(MembraneMeshAttribute::MetaAttribute{s, this}),
    _subSystem(s), _memType(membraneType), _id(_membranes.getID()) {
    
    // Build the meshwork using vertex and triangle information
    _mesh.init(vertexCoordinateList, triangleVertexIndexList);

    size_t numVertices = membraneData.size();
    if(numVertices == 0) return;

    /**************************************************************************
        Setting up vertices and neighbors
    **************************************************************************/
    // Add the vertices
    size_t vertexIndex = 0;
    _vertexVector.reserve(numVertices);
    for(auto& vertexData : membraneData) {
        Vertex* lastAddedVertex = _subSystem->addTrackable<Vertex>(
            array2Vector<double, 3>(get<0>(vertexData)),
            this,
            get<1>(vertexData).size()
        );
        _vertexVector.push_back(lastAddedVertex); // Add to its own storage
        lastAddedVertex->_membraneVertexIdx = vertexIndex++;
    }

    // Register the neighbors
    for(int idx = 0; idx < numVertices; ++idx) {
        auto& neighborData = get<1>(membraneData[idx]);
        Vertex* centerVertex = _vertexVector[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(size_t nIdx = 0; nIdx < numNeighbors; ++nIdx) {
            Vertex* nVertex = _vertexVector[neighborData[nIdx]];
            centerVertex->getNeighborVertices()[nIdx] = nVertex;
            centerVertex->getNeighborVertexIndices()[nVertex] = nIdx;
        }
    }

    /**************************************************************************
        Setting up edges
    **************************************************************************/
    for(int idx = 0; idx < numVertices; ++idx) {
        Vertex* centerVertex = _vertexVector[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(int nIdx = 0; nIdx < numNeighbors; ++nIdx) {

            Vertex* nVertex = centerVertex->getNeighborVertices()[nIdx];

            // Edge registration
            if(centerVertex->getNeighborEdges()[nIdx] == nullptr) { // Edge not registered
                Edge* lastAddedEdge = _subSystem->addTrackable<Edge>(this, centerVertex, nVertex);
                _edgeVector.push_back(lastAddedEdge); // Add to its own storage

                // Bind the edge to vertices, and check whether neighbor exists
                size_t backToCenterIdx = 0;
                try {
                    backToCenterIdx = nVertex->getNeighborVertexIndices().at(centerVertex);
                }
                catch(const std::out_of_range& oor) {
                    cout << "An error occured when trying to add edges of the meshwork. "
                         << "Neighbors must pair with each other. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }

                centerVertex->getNeighborEdges()[nIdx] = lastAddedEdge;
                nVertex->getNeighborEdges()[backToCenterIdx] = lastAddedEdge;
                
                centerVertex->getEdgeHead()[nIdx] = 0;
                nVertex->getEdgeHead()[backToCenterIdx] = 1;

                // Calculate the length of the edge
                lastAddedEdge->getGEdge()->calcLength();
            }
        }
    }

    /**************************************************************************
        Setting up triangles
    **************************************************************************/
    for(int idx = 0; idx < numVertices; ++idx) {
        Vertex* centerVertex = _vertexVector[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(int nIdx = 0; nIdx < numNeighbors; ++nIdx) {

            Vertex* nVertex = centerVertex->getNeighborVertices()[nIdx];
            Vertex* nnVertex = centerVertex->getNeighborVertices()[(nIdx + 1) % numNeighbors];

            // Triangle registration
            if(centerVertex->getNeighborTriangles()[nIdx] == nullptr) { // Triangle not registered
                Triangle* lastAddedTriangle = _subSystem->addTrackable<Triangle>(this, centerVertex, nVertex, nnVertex);
                _triangleVector.push_back(lastAddedTriangle); // Add to its own storage

                // Bind the triangle to vertices, and check whether neighbor exists
                size_t idx12 = 0, idx20 = 0;
                try {
                    idx12 = nVertex->getNeighborVertexIndices().at(nnVertex);
                    idx20 = nnVertex->getNeighborVertexIndices().at(centerVertex);
                }
                catch(const std::out_of_range& oor) {
                    cout << "An error occured when trying to add triangles of the meshwork. "
                         << "Vertices should be able to form triangles. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }

                centerVertex->getNeighborTriangles()[nIdx] = lastAddedTriangle;
                nVertex->getNeighborTriangles()[idx12] = lastAddedTriangle;
                nnVertex->getNeighborTriangles()[idx20] = lastAddedTriangle;
                
                centerVertex->getTriangleHead()[nIdx] = 0;
                nVertex->getTriangleHead()[idx12] = 2;
                nnVertex->getTriangleHead()[idx20] = 1;

                // Bind edges to the triangle
                lastAddedTriangle->getEdges() = {{
                    centerVertex->getNeighborEdges()[nIdx],
                    nVertex->getNeighborEdges()[idx12],
                    nnVertex->getNeighborEdges()[idx20]
                }};
                lastAddedTriangle->getEdgeHead() = {{
                    centerVertex->getEdgeHead()[nIdx],
                    nVertex->getEdgeHead()[idx12],
                    nnVertex->getEdgeHead()[idx20]
                }};
                
                // Bind the triangle to edges
                for(size_t eIdx = 0; eIdx < 3; ++eIdx)
                    lastAddedTriangle->getEdges()[eIdx]->getTriangles()[lastAddedTriangle->getEdgeHead()[eIdx]] = lastAddedTriangle;

                // Calculate the area of the triangle and set it as eqArea
                lastAddedTriangle->getGTriangle()->calcArea();
#ifdef MECHANICS
                lastAddedTriangle->getMTriangle()->setEqArea(
                    lastAddedTriangle->getGTriangle()->getArea() *
                    SysParams::Mechanics().MemEqAreaFactor[membraneType]
                );
#endif
				// Calculate angles for the use of Voronoi cells
				lastAddedTriangle->getGTriangle()->calcTheta();
            }
        }
    }

    /**************************************************************************
        Setting up Voronoi cells
    **************************************************************************/
    for(size_t idx = 0; idx < numVertices; ++idx) {
        GVoronoiCell* gvc = _vertexVector[idx]->getGVoronoiCell();
        gvc->calcArea();
#ifdef MECHANICS
        MVoronoiCell* mvc = _vertexVector[idx]->getMVoronoiCell();
        // Set the current area as eqArea
        mvc->setEqArea(gvc->getArea() * SysParams::Mechanics().MemEqAreaFactor[membraneType]);
#endif
    }

    /**************************************************************************
        Setting up MMembrane object and find volume
    **************************************************************************/
    _gMembrane = unique_ptr<GMembrane>(new GMembrane);
    _gMembrane->setMembrane(this);
    _gMembrane->calcVolume();
#ifdef MECHANICS
    _mMembrane = unique_ptr<MMembrane>(new MMembrane);
    _mMembrane->setMembrane(this);
    _mMembrane->setEqVolume(_gMembrane->getVolume());
#endif

}

void Membrane::printSelf()const {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    cout << "Membrane Id = " << _id << endl;
    cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    cout << "Triangle information..." << endl;
    
    for(auto t : _triangleVector)
        t->printSelf();
    
    cout << endl;
    
}

void Membrane::updateGeometryValue() {
    const auto& vertices = _mesh.getVertices();
    const auto& halfEdges = _mesh.getHalfEdges();
    const auto& edges = _mesh.getEdges();
    const auto& triangles = _mesh.getTriangles();

    const size_t numVertices = vertices.size();
    const size_t numHalfEdges = halfEdges.size();
    const size_t numEdges = edges.size();
    const size_t numTriangles = triangles.size();

    // Calculate angles stored in half edges
    for(size_t hei = 0; hei < numHalfEdges; ++hei) {
        // The angle is (v0, v1, v2)
        const size_t vi0 = _mesh.target(_mesh.prev(hei));
        const size_t vi1 = _mesh.target(hei);
        const size_t vi2 = _mesh.target(_mesh.next(hei));
        const auto& c0 = vertices[vi0].attr.vertex->coordinate;
        const auto& c1 = vertices[vi1].attr.vertex->coordinate;
        const auto& c2 = vertices[vi2].attr.vertex->coordinate;
        auto& heag = _mesh.getHalfEdgeAttribute(hei).gHalfEdge;

        const auto vp = vectorProduct(c1, c0, c1, c2);
        const auto sp = scalarProduct(c1, c0, c1, c2);
        const auto ct = heag.cotTheta = sp / magnitude(vp);
        heag.theta = M_PI_2 - atan(ct);
    }

    // Calculate triangle area, unit normal and cone volume
    for(size_t ti = 0; ti < numTriangles; ++ti) {
        const size_t hei = triangles[ti].halfEdgeIndex;
        const size_t vi0 = _mesh.target(hei);
        const size_t vi1 = _mesh.target(_mesh.next(hei));
        const size_t vi2 = _mesh.target(_mesh.prev(hei));
        const auto& c0 = vertices[vi0].attr.vertex->coordinate;
        const auto& c1 = vertices[vi1].attr.vertex->coordinate;
        const auto& c2 = vertices[vi2].attr.vertex->coordinate;
        auto& tag = _mesh.getTriangleAttribute(ti).gTriangle;

        const auto vp = vectorProduct(c0, c1, c0, c2);

        // area
        tag.area = magnitude(vp) * 0.5;

        // unit normal
        tag.unitNormal = vector2Vec<3, double>(normalize(vp));

        // cone volume
        tag.coneVolume = dotProduct(c0, vp) / 6;
    }

    // Calculate edge length and pesudo unit normal
    for(size_t ei = 0; ei < numEdges; ++ei) {
        const size_t hei = edges[ei].halfEdgeIndex;
        const size_t vi0 = _mesh.target(hei);
        const size_t vi1 = _mesh.target(_mesh.prev(hei));

        // length
        _mesh.getEdgeAttribute(ei).gEdge.length = twoPointDistance(vertices[vi0].attr.vertex->coordinate, vertices[vi1].attr.vertex->coordinate);

        // pseudo unit normal
        if(halfEdges[hei].hasOpposite) {
            const size_t ti0 = _mesh.triangle(hei);
            const size_t ti1 = _mesh.triangle(_mesh.opposite(hei));
            _mesh.getEdgeAttribute(ei).gEdge.pseudoUnitNormal = normalize(
                triangles[ti0].attr.gTriangle.unitNormal + triangles[ti1].attr.gTriangle.unitNormal
            );
        }
    }

    // Calculate vcell area, curvature and vertex pseudo unit normal
    for(size_t vi = 0; vi < numVertices; ++vi) {
        auto& vag = _mesh.getVertexAttribute(vi).gVertex;

        // clearing
        vag.area = 0.0;
        vag.pseudoUnitNormal = {0.0, 0.0, 0.0};

        // k1 = 2A * k, where k is the result of LB operator
        Vec3 k1 {};

        _mesh.forEachHalfEdgeTargetingVertex(vi, [&mesh = _mesh, &vertices, &vi, &vag, &k1](size_t hei) {
            const size_t hei_o = mesh.opposite(hei);
            const size_t ti0 = mesh.triangle(hei);
            const size_t ti1 = mesh.triangle(hei_o);
            const size_t vn = mesh.target(hei_o);
            const size_t hei_n = mesh.next(hei);
            const size_t hei_on = mesh.next(hei_o);
            const auto ci = vector2Vec<3, double>(vertices[vi].attr.vertex->coordinate);
            const auto cn = vector2Vec<3, double>(vertices[vn].attr.vertex->coordinate);

            const auto sumCotTheta = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.cotTheta + mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.cotTheta;

            const auto theta = mesh.getHalfEdgeAttribute(hei).gHalfEdge.theta;

            const Vec3 diff = ci - cn;
            const auto dist2 = magnitude2(diff);

            vag.area += sumCotTheta * dist2 * 0.125;

            k1 += sumCotTheta * diff;
            vag.pseudoUnitNormal += theta * mesh.getTriangleAttribute(ti0).gTriangle.unitNormal;
        });

        const double invA = 1 / vag.area;
        const double magK1 = magnitude(k1);

        normalize(vag.pseudoUnitNormal);

        const int flippingCurv = (dotProduct(k1, vag.pseudoUnitNormal) > 0 ? 1 : -1);

        vag.curv = fliipingCurv * magK1 * 0.25 * invA;
    }
}
void Membrane::updateGeometryValueWithDerivative() {

    const auto& vertices = _mesh.getVertices();
    const auto& halfEdges = _mesh.getHalfEdges();
    const auto& edges = _mesh.getEdges();
    const auto& triangles = _mesh.getTriangles();

    const size_t numVertices = vertices.size();
    const size_t numHalfEdges = halfEdges.size();
    const size_t numEdges = edges.size();
    const size_t numTriangles = triangles.size();

    // Calculate edge length with deriviative
    for(size_t ei = 0; ei < numEdges; ++ei) {
        // The edge must have an opposite
        const size_t hei = edges[ei].halfEdgeIndex;
        const size_t hei_o = _mesh.opposite(hei);
        const size_t vi0 = _mesh.target(hei);
        const size_t vi1 = _mesh.target(hei_o);
        const auto c0 = vector2Vec<3, double>(vertices[vi0].attr.vertex->coordinate);
        const auto c1 = vector2Vec<3, double>(vertices[vi1].attr.vertex->coordinate);

        const auto length = _mesh.getEdgeAttribute(ei).gEdge.length = dist(c0, c1);
        _mesh.getHalfEdgeAttribute(hei).gHalfEdge.dEdgeLength = (c0 - c1) / length;
        _mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dEdgeLength = (c1 - c0) / length;
    }

    // Calculate angles and triangle areas with derivative
    // Calculate triangle normals and cone volumes
    for(size_t ti = 0; ti < numTriangles; ++ti) {
        size_t hei[3];
        hei[0] = triangles[ti].halfEdgeIndex;
        hei[1] = _mesh.next(hei[0]);
        hei[2] = _mesh.next(hei[1]);
        auto& tag = _mesh.getTriangleAttribute(ti).gTriangle;

        const size_t vi[] {_mesh.target(hei[0]), _mesh.target(hei[1]), _mesh.target(hei[2])};
        const Vec3 c[] {
            vector2Vec<3, double>(vertices[vi[0]].attr.vertex->coordinate),
            vector2Vec<3, double>(vertices[vi[1]].attr.vertex->coordinate),
            vector2Vec<3, double>(vertices[vi[2]].attr.vertex->coordinate)
        };

        const double l[] {
            edges[_mesh.edge(hei[0])].attr.length,
            edges[_mesh.edge(hei[1])].attr.length,
            edges[_mesh.edge(hei[2])].attr.length
        };

        const double dots[] {
            dot(c[1] - c[0], c[2] - c[0]),
            dot(c[2] - c[1], c[0] - c[1]),
            dot(c[0] - c[2], c[1] - c[2])
        };

        const auto cp = cross(c[1] - c[0], c[2] - c[0]); // Pointing outward

        const auto r0 = cross(c[1], c[2]); // Used in cone volume

        // Calculate area
        const auto area = tag.area = magnitude(cp) * 0.5;
        const auto invA = 1.0 / area;

        // Calculate area gradients
        {
            const auto r01 = c[1] - c[0];
            const auto r02 = c[2] - c[0];
            _mesh.getHalfEdgeAttribute(hei[0]).gHalfEdge.dArea = (-l[1]*l[1]* r02 - l[0]*l[0]* r01 + dots[0]*(r01 + r02)) * (invA * 0.25);
            _mesh.getHalfEdgeAttribute(hei[1]).gHalfEdge.dArea = (l[0]*l[0]* r01 - dots[0]* r02) * (invA * 0.25);
            _mesh.getHalfEdgeAttribute(hei[2]).gHalfEdge.dArea = (l[1]*l[1]* r02 - dots[0]* r01) * (invA * 0.25);
        }

        // Calculate thetas and gradients
        for(size_t ai = 0; ai < 3; ++ai) {
            auto& heag = _mesh.getHalfEdgeAttribute(hei[ai]).gHalfEdge;

            const auto ct = heag.cotTheta = dots[ai] * invA * 0.5;
            heag.theta = M_PI_2 - atan(ct);

            const size_t ai_n = (ai + 1) % 3;
            const sizs_t ai_p = (ai + 2) % 3;
            auto& heag_n = _mesh.getHalfEdgeAttribute(hei[ai_n]).gHalfEdge;
            auto& heag_p = _mesh.getHalfEdgeAttribute(hei[ai_p]).gHalfEdge;

            const auto r01 = c[ai_n] - c[ai];
            const auto r02 = c[ai_p] - c[ai];

            heag.dCotTheta[1] =
                -(r01 + r02) * (invA * 0.5)
                -(dots[ai] * invA * invA * 0.5) * heag.dArea;
            heag.dCotTheta[2] = r02 * (invA * 0.5) - (dots[ai] * invA * invA * 0.5) * heag_n.dArea;
            heag.dCotTheta[0] = r01 * (invA * 0.5) - (dots[ai] * invA * invA * 0.5) * heag_p.dArea;
        }

        // Calculate unit normal
        tag.unitNormal = normalizedVector(cp);

        // Calculate cone volume and derivative
        tag.coneVolume = dot(c[0], r0) / 6.0;
        // The derivative of cone volume will be accumulated to each vertex

    }

    // Clear derivative of vcell area on neighbors
    for(size_t hei = 0; hei < numHalfEdges; ++hei) {
        _mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborArea = {0.0, 0.0, 0.0};
    }

    // Calculate vcell area, curvature with derivative
    // Calculate vertex pseudo unit normal
    // Calculate derivative of volume on vertices
    for(size_t vi = 0; vi < numVertices; ++vi) {
        auto& vag = _mesh.getVertexAttribute(vi).gVertex;

        // clearing
        vag.area = 0.0;
        vag.dArea = {0.0, 0.0, 0.0};
        vag.pseudoUnitNormal = {0.0, 0.0, 0.0};
        vag.dVolume = {0.0, 0.0, 0.0};

        // K = 2*H*n is the result of LB operator
        // And let k1 = 2*A*k (as an intermediate variable)
        Vec3 k1 {};

        // derivative of k1 and curvature will be calculated in the next loop
        _mesh.forEachHalfEdgeTargetingVertex(vi, [&mesh = _mesh, &vertices, &vi, &vag, &k1](size_t hei) {
            const size_t hei_o = mesh.opposite(hei);
            const size_t ti0 = mesh.triangle(hei);
            const size_t ti1 = mesh.triangle(hei_o);
            const size_t vn = mesh.target(hei_o);
            const size_t hei_n = mesh.next(hei);
            const size_t hei_on = mesh.next(hei_o);
            const size_t hei_right = mesh.opposite(mesh.next(hei_on)); // hei_left is hei_n
            const auto ci = vector2Vec<3, double>(vertices[vi].attr.vertex->coordinate);
            const auto cn = vector2Vec<3, double>(vertices[vn].attr.vertex->coordinate);

            const auto sumCotTheta = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.cotTheta + mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.cotTheta;
            const auto& dCotThetaLeft = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.dCotTheta;
            const auto& dCotThetaRight = mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.dCotTheta;

            const auto theta = mesh.getHalfEdgeAttribute(hei).gHalfEdge.theta;

            const Vec3 diff = ci - cn;
            const auto dist2 = magnitude2(diff);

            vag.area += sumCotTheta * dist2 * 0.125;

            // Area derivative
            vag.dArea +=
                (dCotThetaLeft[0] + dCotThetaRight[2]) * (dist2 * 0.125)
                + (sumCotTheta * 0.25) * diff; // d(dist2) / dx = 2 * diff
            mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborArea +=
                (dCotThetaLeft[2] + dCotThetaRight[0]) * (dist2 * 0.125)
                - (sumCotTheta * 0.25) * diff; // d(dist2) / dx = -2 * diff
            mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.dNeighborArea +=
                dCotThetaLeft[1] * (dist2 * 0.125);
            mesh.getHalfEdgeAttribute(hei_right).gHalfEdge.dNeighborArea +=
                dCotThetaRight[1] * (dist2 * 0.125);

            k1 += sumCotTheta * diff;
            vag.pseudoUnitNormal += theta * mesh.getTriangleAttribute(ti0).gTriangle.unitNormal;
        });

        // Calculate derivative of k1 and curvature
        // Using another loop because k1 is needed for curvature derivative
        Mat3 dK1 {}; // On center vertex
        _mesh.forEachHalfEdgeTargetingVertex(vi, [this, vi, &vag, &dK1](size_t hei) {
            const size_t hei_o = mesh.opposite(hei);
            const size_t vn = mesh.target(hei_o);
            const size_t hei_n = mesh.next(hei);
            const size_t hei_on = mesh.next(hei_o);
            const size_t hei_right = mesh.opposite(mesh.next(hei_on)); // hei_left is hei_n
            const auto ci = vector2Vec<3, double>(vertices[vi].attr.vertex->coordinate);
            const auto cn = vector2Vec<3, double>(vertices[vn].attr.vertex->coordinate);

            // Accumulate dK1 on the center vertex
            dK1 += tensorTmp0 + Eye3CArray * sumCotTheta;
            where tensorTmp0 = (<sum-of-d-cot-theta>)T * diff;

            // Calculate dK1 on target vertex
            // As direct target
        });
        // TODO

            // Temporary variables
            double matTmp1[9], matTmp2[9];
            double vecTmp[3];

            // NOTE: Not necessarily correct for peripheral vertices
            for(size_t nIdx = 0; nIdx < nN; ++nIdx) {

                double diff[3];

                // Calculate area, K1 (= 2A * 2H * normal) and their derivatives. Also calcs pseudo unit normal w/o normalization
                arrArea[idx] += sumCotTheta * dist2 * 0.125;

                // Now, dDiff is Eye3, and dNDiff is -Eye3
                double tensorTmp0[9];
                double tensorTmp1[9];
                double tensorTmpL[9];
                double tensorTmpR[9];
                tensorProduct(tensorTmp0, vectorSum(vecTmp, arrDCotThetaTriangle + 27*triLIdx + 9*triLIdx1 + 3*triLIdx0, arrDCotThetaTriangle + 27*triRIdx + 9*triRIdx2 + 3*triRIdx0), diff);
                tensorProduct(tensorTmp1, vectorSum(vecTmp, arrDCotThetaTriangle + 27*triLIdx + 9*triLIdx1 + 3*triLIdx2, arrDCotThetaTriangle + 27*triRIdx + 9*triRIdx2 + 3*triRIdx1), diff);
                tensorProduct(tensorTmpL, arrDCotThetaTriangle + 27*triLIdx + 9*triLIdx1 + 3*triLIdx1, diff);
                tensorProduct(tensorTmpR, arrDCotThetaTriangle + 27*triRIdx + 9*triRIdx2 + 3*triRIdx2, diff);

                matrixIncrease(dK1, matrixSum(matTmp2, tensorTmp0, matrixMultiply(matTmp1, Eye3CArray, sumCotTheta)));
                matrixIncrease(dNeighborK1.get() + 9*nIdx, matrixDifference(matTmp2, tensorTmp1, matrixMultiply(matTmp1, Eye3CArray, sumCotTheta)));
                matrixIncrease(dNeighborK1.get() + 9*nIdxP, tensorTmpL);
                matrixIncrease(dNeighborK1.get() + 9*nIdxN, tensorTmpR);

                // Added to derivative of sum of cone volume
                crossProduct(vecTmp, arrCoord + 3*arrVCellTopoInfo[idx0BIdx + 1 + nIdx], arrCoord + 3*arrVCellTopoInfo[idx0BIdx + 1 + nIdxN]);
                arrDVolumeVertex[3*idx + 0] += vecTmp[0] / 6;
                arrDVolumeVertex[3*idx + 1] += vecTmp[1] / 6;
                arrDVolumeVertex[3*idx + 2] += vecTmp[2] / 6;
            } // End loop neighbors

            const double invA = 1 / arrArea[idx];
            const double magK1 = magnitude(k1);

            // Calculate pseudo unit normal
            normalize(arrPseudoUnitNormal + 3*idx);

            // Calculate mean curvature H = |k1| / 4A
            // dH = (dK1)K1 / 4A|K1| - |K1|dA / 4A^2
            const int flippingCurv = (dotProduct(k1, arrPseudoUnitNormal + 3*idx) > 0 ? 1 : -1);
            const double dCurvFac1 = 0.25 * invA / magK1;
            const double dCurvFac2 = -0.25 * invA * invA * magK1;

            matrixVectorProduct(vecTmp, dK1, k1);
            for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                arrDCurv[3*bsIdx + coordIdx] = flippingCurv * (dCurvFac1 * vecTmp[coordIdx] + dCurvFac2 * arrDArea[3*bsIdx + coordIdx]);
            }
            for(size_t nIdx = 0; nIdx < nN; ++nIdx) {
                matrixVectorProduct(vecTmp, dNeighborK1.get() + 9*nIdx, k1);
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    arrDCurv[3*bsIdx + 3*(nIdx+1) + coordIdx] = flippingCurv * (dCurvFac1 * vecTmp[coordIdx] + dCurvFac2 * arrDArea[3*bsIdx + 3*(nIdx+1) + coordIdx]);
                }
            }
            arrCurv[idx] = flippingCurv * magK1 * 0.25 * invA;

    } // End loop vertices (V cells)
}

void Membrane::updateGeometry(bool calcDerivative, double d) {
    /**************************************************************************
    Updates the geometric properties of all the elements of the membrane. This
    MUST be called before any energy calculation of the membrane.

    Due to the geometry dependencies, the operation order below is important.

    Some of the geometric properties are not needed in energy computation, so
    for efficiency, they are not updated in this function. They must be updated
    manually before use. They are:
        - <None>
    **************************************************************************/

    for(auto& e: _edgeVector) {
        auto ge = e->getGEdge();
        if(calcDerivative) ge->calcLength(); else ge->calcStretchedLength(d);
    }
    for(auto& t: _triangleVector) {
        auto gt = t->getGTriangle();
        if(calcDerivative) gt->calcTheta(); else gt->calcStretchedTheta(d);
        if(calcDerivative) gt->calcArea(); else gt->calcStretchedArea(d);
        if(calcDerivative) gt->calcUnitNormal(); else gt->calcStretchedUnitNormal(d);
    }
    for(auto& e: _edgeVector) {
        auto ge = e->getGEdge();
        if(calcDerivative) ge->calcPseudoUnitNormal(); else ge->calcStretchedPseudoUnitNormal(d);
    }
    for(auto& v: _vertexVector) {
        auto gv = v->getGVoronoiCell();
        if(calcDerivative) gv->calcArea(); else gv->calcStretchedArea(d);
        if(calcDerivative) gv->calcCurv(); else gv->calcStretchedCurv(d);
        if(calcDerivative) gv->calcPseudoUnitNormal(); else gv->calcStretchedPseudoUnitNormal(d);
    }
    if(calcDerivative) _gMembrane->calcVolume(); else _gMembrane->calcStretchedVolume(d);
}

double Membrane::signedDistance(const std::array<double, 3>& p, bool safe)const {
    if(!_isClosed) throw std::logic_error("Membrane is not closed while trying to find signed distance field.");

    /**************************************************************************
    The function works in the following procedure:

    - Iterate through a certain amount of triangles, and for each triangle
        - Find the projection of the point on the triangle plane, and determine
          which element is responsible for being the closest to the point
          (which vertex/ which edge/ this triangle).
        - Find the unsigned distance with the closest element, and then find
          the signed distance using the normal or pseudo normal. Record the
          value with the smallest unsigned distance.
    
    Before this function is used, the following must be calculated:
        - The positions of all the elements are updated
        - The normal and pseudo normal at the triangles, edges and vertices
        - The length of edges
    
    In fact, the signed distance field serves as a good candidate for membrane
    boundary potential. However, this field is only C0-continuous, so some
    configurations might result in non-stop situation in CG method.
    **************************************************************************/

    /**************************************************************************
    Determine the triangles needed to be looped through

    Safe        : Use all triangles
    Not safe    : Search neighboring 27 compartments for triangles
    **************************************************************************/
    const vector<Triangle*>* loopingTriangles = nullptr;

    unique_ptr<vector<Triangle*>> limitedLoopingTriangles;
    if(safe) {
        // Loop through all the triangles
        loopingTriangles = &_triangleVector;
    }
    else {
        // Find the current compartment indices containing point p
        vector<size_t> indices = GController::getCompartmentIndices(mathfunc::array2Vector(p));

        if(indices.size() != 3)
            throw std::logic_error("Only 3D compartments are allowed in the membrane signed distance calculation.");

		unordered_set<Triangle*> triSet; // Set of triangles to be searched

        // Find all the neighbor compartments and find triangles to be added to search set.
        // The for loops are written like this for aesthetic reasons. Be careful with the scope below.
        for(int v0 = -1; v0 <= 1; ++v0) for(int v1 = -1; v1 <= 1; ++v1) for(int v2 = -1; v2 <= 1; ++v2) {

            vector<size_t> newIndices = {indices[0] + v0, indices[1] + v1, indices[2] + v2};
            if(GController::indicesOutOfBound(newIndices)) {
                // Compartment not found. Simply ignore it.
                continue;
            }

            Compartment* c = GController::getCompartment(newIndices);

            // Compartment exists
            const unordered_set<Triangle*>& triSetEach = c->getTriangles();
            triSet.insert(triSetEach.begin(), triSetEach.end());
        }

        if(triSet.empty()) {
            cout << "Warning: triangles not found in neighboring compartments. "
                 << "Getting the result from safe mode."
                 << endl;
            return signedDistance(p, true);
        }

        // Copy the triangles in the unordered set into a vector
        limitedLoopingTriangles = unique_ptr<vector<Triangle*>>(new vector<Triangle*>(triSet.begin(), triSet.end()));

        // Make the new vector the vector to be looped
        loopingTriangles = limitedLoopingTriangles.get();
    }

    double minAbsDistance = numeric_limits<double>::infinity();
    for(Triangle* t: *loopingTriangles) {
        /**********************************************************************
        Calculate the barycentric coordinate of the projection point p'

        See Heidrich 2005, Computing the Barycentric Coordinates of a Projected
        Point.
        **********************************************************************/
        const array<double, 3> v0 = vector2Array<double, 3>(t->getVertices()[0]->coordinate);
        const array<double, 3> v1 = vector2Array<double, 3>(t->getVertices()[1]->coordinate);
        const array<double, 3> v2 = vector2Array<double, 3>(t->getVertices()[2]->coordinate);

        array<double, 3> n = vectorProduct(v0, v1, v0, v2);
        double oneOver4AreaSquared = 1.0 / dotProduct(n, n);

        double b1 = dotProduct(vectorProduct(v0, p, v0, v2), n) * oneOver4AreaSquared;
        double b2 = dotProduct(vectorProduct(v0, v1, v0, p), n) * oneOver4AreaSquared;
        double b0 = 1.0 - b1 - b2;

        // Now p' = b0*v0 + b1*v1 + b2*v2

        double d = numeric_limits<double>::infinity();
        if(b0 >= 0 && b1 >= 0 && b2 >= 0) {
            // p' is inside the triangle
            d = dotProduct(t->getGTriangle()->getUnitNormal(), vectorDifference(p, v0));
        } else {
            // p' is outside the triangle
            const array<double, 3> r {
                t->getEdges()[1]->getGEdge()->getLength(), // 1->2
                t->getEdges()[2]->getGEdge()->getLength(), // 2->0
                t->getEdges()[0]->getGEdge()->getLength()  // 0->1
            };
            double dot_1p_12 = scalarProduct(v1, p, v1, v2);
            double dot_2p_20 = scalarProduct(v2, p, v2, v0);
            double dot_0p_01 = scalarProduct(v0, p, v0, v1);

            if(b0 < 0 && dot_1p_12 >= 0 && dot_1p_12 <= r[0]*r[0]) {
                // On edge 12
                d = magnitude(vectorProduct(v1, p, v1, v2)) / r[0];
                if(dotProduct(t->getEdges()[1]->getGEdge()->getPseudoUnitNormal(), vectorDifference(p, v1)) < 0) d = -d;
            } else if(b1 < 0 && dot_2p_20 >= 0 && dot_2p_20 <= r[1]*r[1]) {
                // On edge 20
                d = magnitude(vectorProduct(v2, p, v2, v0)) / r[1];
                if(dotProduct(t->getEdges()[2]->getGEdge()->getPseudoUnitNormal(), vectorDifference(p, v2)) < 0) d = -d;
            } else if(b2 < 0 && dot_0p_01 >= 0 && dot_0p_01 <= r[2]*r[2]) {
                // On edge 01
                d = magnitude(vectorProduct(v0, p, v0, v1)) / r[2];
                if(dotProduct(t->getEdges()[0]->getGEdge()->getPseudoUnitNormal(), vectorDifference(p, v0)) < 0) d = -d;
            } else if(dot_0p_01 < 0 && dot_2p_20 > r[1]*r[1]) {
                // On vertex 0
                d = twoPointDistance(v0, p);
                if(dotProduct(t->getVertices()[0]->getGVoronoiCell()->getPseudoUnitNormal(), vectorDifference(p, v0)) < 0)
                    d = -d;
            } else if(dot_1p_12 < 0 && dot_0p_01 > r[2]*r[2]) {
                // On vertex 1
                d = twoPointDistance(v1, p);
                if(dotProduct(t->getVertices()[1]->getGVoronoiCell()->getPseudoUnitNormal(), vectorDifference(p, v1)) < 0)
                    d = -d;
            } else if(dot_2p_20 < 0 && dot_1p_12 > r[0]*r[0]) {
                // On vertex 2
                d = twoPointDistance(v2, p);
                if(dotProduct(t->getVertices()[2]->getGVoronoiCell()->getPseudoUnitNormal(), vectorDifference(p, v2)) < 0)
                    d = -d;
            } else {
                // The program should never go here
                throw logic_error("Unknown case of point projection on the plane of triangle.");
            }
        }

        // Update with distance with less absolute value
        if(abs(d) < abs(minAbsDistance)) minAbsDistance = d;
    }
    
    return minAbsDistance;
}

bool Membrane::contains(const std::array<double, 3>& point)const {
    return signedDistance(point, true) < 0.0;
}

double Membrane::meshworkQuality()const {
    /*
    This function calculates the quality of the meshwork of this membrane, and
    the result is represented as a value between 0 and 1, 0 being worst and 1
    being the best case (all equilateral).

    The criterion used is the minimum angle in each triangle, parametrized and
    averaged over all triangles. And the calculation requires the result of
        - The angle calculation of all triangles
    */

    double res = 0;

    for(Triangle* t: _triangleVector) {
        auto gt = t->getGTriangle();

        // Find minimum angle
        double minAngle = M_PI / 3; // The largest possible minimum angle in a triangle
        for(double eachTheta: gt->getTheta())
            if(eachTheta < minAngle) minAngle = eachTheta;
        
        // Parametrize the minimum angle and add to result
        double q = minAngle * 3 / M_PI; // Normalize, so that q should be within 0 and 1
        q *= q; // Small angles get more "emphasized"
        res += q;
    }

    res /= _triangleVector.size(); // Average instead of sum

    return res;
}
