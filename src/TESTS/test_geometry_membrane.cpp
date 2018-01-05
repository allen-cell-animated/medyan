
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifdef TESTING

#  define DO_THIS_GEOMETRY_MEMBRANE_TEST
#  ifdef DO_THIS_GEOMETRY_MEMBRANE_TEST

#    include "gtest/gtest.h"

#    include <random>

#    include "common.h"
#    include "MathFunctions.h"
using namespace mathfunc;
#    include "Rand.h"

#    include "GController.h"
#    include "SubSystem.h"

#    include "Vertex.h"
#    include "Edge.h"
#    include "Triangle.h"
#    include "MVoronoiCell.h"
#    include "GEdge.h"
#    include "GTriangle.h"
#    include "Membrane.h"

namespace {
    using VertexData = tuple<array<double, 3>, vector<size_t>>;
    using MembraneData = vector<VertexData>;

    MembraneData membraneDataOctahedron(double radius) {
        double c = 2 * radius;

        return MembraneData{
            VertexData({c, c, c+radius}, {1, 2, 3, 4}),
            VertexData({c+radius, c, c}, {0, 4, 5, 2}),
            VertexData({c, c+radius, c}, {0, 1, 5, 3}),
            VertexData({c-radius, c, c}, {0, 2, 5, 4}),
            VertexData({c, c-radius, c}, {0, 3, 5, 1}),
            VertexData({c, c, c-radius}, {4, 3, 2, 1})
        };
    }

    class MembraneGeometryTest: public ::testing::Test {
    protected:
        double radius;
        SubSystem s;
        MembraneData memData;
        Membrane *m;

        MembraneGeometryTest(): radius(100), memData(membraneDataOctahedron(radius)) {
            SysParams::GParams.compartmentSizeX = 1e10;
            SysParams::GParams.compartmentSizeY = 1e10;
            SysParams::GParams.compartmentSizeZ = 1e10;
            
            SysParams::GParams.NX = 1;
            SysParams::GParams.NY = 1;
            SysParams::GParams.NZ = 1;
            
            SysParams::GParams.nDim = 3;

            GController g(&s); // Dummy variable to initialize the compartments
            g.initializeGrid();

            SysParams::GParams.cylinderNumMon.resize(1, 3);
            m = new Membrane(&s, 0, memData);
        }
        ~MembraneGeometryTest() {
            SysParams::GParams.cylinderNumMon.resize(0);
            delete m;
        }

    };

    void recordCoordinate(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->coordinateP = it->coordinate;
    }
    void resetCoordinate(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->coordinate = it->coordinateP;
    }
    void assignRandomForce(Membrane* m, double sigma) {
        normal_distribution<> nd(0, sigma);

        for(Vertex* it: m->getVertexVector()) {
            for(double& eachForce: it->force) {
                eachForce = nd(Rand::engFixed);
            }
        }
    }
    void moveAlongForce(Membrane* m, double d) {
        for(Vertex* it: m->getVertexVector()) {
            for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                it->coordinate[coordIdx] += it->force[coordIdx] * d;
            }
        }
    }
}

TEST_F(MembraneGeometryTest, Topology) {

    // Check that the vertices, edges and triangles are correctly registered.
    EXPECT_EQ(m->getVertexVector().size(), 6);
    EXPECT_EQ(m->getEdgeVector().size(), 12);
    EXPECT_EQ(m->getTriangleVector().size(), 8);
    
}

TEST_F(MembraneGeometryTest, Geometry) {

    /**************************************************************************
        Check normal geometry
    **************************************************************************/
    m->updateGeometry(true);

    // Check edge length
    double exEdgeLen = radius * sqrt(2);
    for(Edge* it: m->getEdgeVector())
        EXPECT_DOUBLE_EQ(it->getGEdge()->getLength(), exEdgeLen);
    // Check triangle area and angle
    double exTriangleArea = radius * radius * sqrt(3) / 2;
    double exTriangleAngle = M_PI / 3;
    for(Triangle* it: m->getTriangleVector()) {
        EXPECT_DOUBLE_EQ(it->getGTriangle()->getArea(), exTriangleArea);
        for(double eachTheta: it->getGTriangle()->getTheta())
            EXPECT_DOUBLE_EQ(eachTheta, exTriangleAngle);
    }
    // Check Voronoi cell area and curvature
    double exVCellArea = radius * radius * sqrt(3) * 2 / 3;
    double exVCellCurv = 1 / radius;
    for(Vertex* it: m->getVertexVector()) {
        EXPECT_DOUBLE_EQ(it->getMVoronoiCell()->getArea(), exVCellArea);
        EXPECT_DOUBLE_EQ(it->getMVoronoiCell()->getCurv(), exVCellCurv);
    }

    /**************************************************************************
        Check stretched geometry
    **************************************************************************/
    Vertex* v0 = m->getVertexVector()[0];
    size_t numNeighborV0 = v0->getNeighborNum();
    ASSERT_EQ(numNeighborV0, 4) << "The 0th vertex does not have the right number of neighbors.";

    // Assign force to 0th vertex
    v0->force[2] = -radius;

    m->updateGeometry(false, 1.0); // Moves the 0th vertex to the center of the octahedron.

    // Check edge length
    double exStretchedEdgeLen = radius;
    for(Edge* it: v0->getNeighborEdges())
        EXPECT_DOUBLE_EQ(it->getGEdge()->getStretchedLength(), exStretchedEdgeLen);
    // Check triangle area and angle
    double exStretchedTriangleArea = radius * radius / 2;
    double exStretchedTriangleAngleIn = M_PI / 2;
    double exStretchedTriangleAngleOut = M_PI / 4;
    for(size_t nIdx = 0; nIdx < numNeighborV0; ++nIdx) {
        EXPECT_DOUBLE_EQ(v0->getNeighborTriangles()[nIdx]->getGTriangle()->getStretchedArea(), exStretchedTriangleArea);
        for(size_t angleIdx = 0; angleIdx < 3; ++angleIdx) {
            if((3 - v0->getTriangleHead()[nIdx]) % 3 == angleIdx)
                EXPECT_DOUBLE_EQ(
                    v0->getNeighborTriangles()[nIdx]->getGTriangle()->getStretchedTheta()[angleIdx],
                    exStretchedTriangleAngleIn
                );
            else
                EXPECT_DOUBLE_EQ(
                    v0->getNeighborTriangles()[nIdx]->getGTriangle()->getStretchedTheta()[angleIdx],
                    exStretchedTriangleAngleOut
                );
        }
    }
    // Check Voronoi cell area and curvature
    double exStretchedVCellArea = radius * radius;
    double exStretchedVCellCurv = 0;
    EXPECT_DOUBLE_EQ(v0->getMVoronoiCell()->getStretchedArea(), exStretchedVCellArea);
    EXPECT_DOUBLE_EQ(v0->getMVoronoiCell()->getStretchedCurv(), exStretchedVCellCurv);

}

TEST_F(MembraneGeometryTest, Derivative) {
    m->updateGeometry(true);
    recordCoordinate(m);
    assignRandomForce(m, radius/200); // Simple test shows that 100 induces a change not small enough

    size_t numEdges = m->getEdgeVector().size();
    size_t numTriangles = m->getTriangleVector().size();
    size_t numVertices = m->getVertexVector().size();

    // Then move every vertex a little bit
    moveAlongForce(m, 1.0);
    m->updateGeometry(false, 0.0); // use stretched calculation here to speed things up
    
    vector<double> edgeLength1(numEdges);
    for(size_t idx = 0; idx < numEdges; ++idx) {
        edgeLength1[idx] = m->getEdgeVector()[idx]->getGEdge()->getStretchedLength();
    }
    vector<double> triangleArea1(numTriangles);
    vector<array<double, 3>> triangleTheta1(numTriangles, {{}});
    vector<array<double, 3>> triangleCotTheta1(numTriangles, {{}});
    for(size_t idx = 0; idx < numTriangles; ++idx) {
        triangleArea1[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedArea();
        triangleTheta1[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedTheta();
        triangleCotTheta1[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedCotTheta();
    }
    vector<double> vCellArea1(numVertices);
    vector<double> vCellCurv1(numVertices);
    for(size_t idx = 0; idx < numVertices; ++idx) {
        vCellArea1[idx] = m->getVertexVector()[idx]->getMVoronoiCell()->getStretchedArea();
        vCellCurv1[idx] = m->getVertexVector()[idx]->getMVoronoiCell()->getStretchedCurv();
    }

    // Now move every vertex in the opposite direction
    resetCoordinate(m);
    moveAlongForce(m, -1.0);
    m->updateGeometry(false, 0.0);

    vector<double> edgeLength2(numEdges);
    for(size_t idx = 0; idx < numEdges; ++idx) {
        edgeLength2[idx] = m->getEdgeVector()[idx]->getGEdge()->getStretchedLength();
    }
    vector<double> triangleArea2(numTriangles);
    vector<array<double, 3>> triangleTheta2(numTriangles, {{}});
    vector<array<double, 3>> triangleCotTheta2(numTriangles, {{}});
    for(size_t idx = 0; idx < numTriangles; ++idx) {
        triangleArea2[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedArea();
        triangleTheta2[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedTheta();
        triangleCotTheta2[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedCotTheta();
    }
    vector<double> vCellArea2(numVertices);
    vector<double> vCellCurv2(numVertices);
    for(size_t idx = 0; idx < numVertices; ++idx) {
        vCellArea2[idx] = m->getVertexVector()[idx]->getMVoronoiCell()->getStretchedArea();
        vCellCurv2[idx] = m->getVertexVector()[idx]->getMVoronoiCell()->getStretchedCurv();
    }

    // Compare the results with derivative predictions
    // A(x+h) - A(x-h) = dotProduct(2h, dA/dx)
    for(size_t idx = 0; idx < numEdges; ++idx) {
        Edge* e = m->getEdgeVector()[idx];
        // Edge length
        double exDiff = 0.0;
        for(size_t vIdx = 0; vIdx < 2; ++vIdx) {
            exDiff += 2 * dotProduct(
                e->getVertices()[vIdx]->force,
                array2Vector<double, 3>(e->getGEdge()->getDLength()[vIdx])
            );
        }
        EXPECT_NEAR(edgeLength1[idx] - edgeLength2[idx], exDiff, abs(exDiff / 1000));
    }
	for(size_t idx = 0; idx < numTriangles; ++idx) {
        Triangle *t = m->getTriangleVector()[idx];
        // Triangle area
        double exDiff = 0.0;
        for(size_t vIdx = 0; vIdx < 3; ++vIdx) {
            exDiff += 2 * dotProduct(
                t->getVertices()[vIdx]->force,
                array2Vector<double, 3>(t->getGTriangle()->getDArea()[vIdx])
            );
        }
        EXPECT_NEAR(triangleArea1[idx] - triangleArea2[idx], exDiff, abs(exDiff / 1000));
        // Triangle angles
        for(size_t aIdx = 0; aIdx < 3; ++aIdx) {
            exDiff = 0.0;
            for(size_t vIdx = 0; vIdx < 3; ++vIdx) {
                exDiff += 2 * dotProduct(
                    t->getVertices()[vIdx]->force,
                    array2Vector<double, 3>(t->getGTriangle()->getDTheta()[aIdx][vIdx])
                );
            }
            EXPECT_NEAR(triangleTheta1[idx][aIdx] - triangleTheta2[idx][aIdx], exDiff, abs(exDiff / 1000));
        }
        // Triangle angles (cot)
        for(size_t aIdx = 0; aIdx < 3; ++aIdx) {
            exDiff = 0.0;
            for(size_t vIdx = 0; vIdx < 3; ++vIdx) {
                exDiff += 2 * dotProduct(
                    t->getVertices()[vIdx]->force,
                    array2Vector<double, 3>(t->getGTriangle()->getDCotTheta()[aIdx][vIdx])
                );
            }
            EXPECT_NEAR(triangleCotTheta1[idx][aIdx] - triangleCotTheta2[idx][aIdx], exDiff, abs(exDiff / 1000));
        }
    }
    for(size_t idx = 0; idx < numVertices; ++idx) {
        Vertex *v = m->getVertexVector()[idx];
        size_t numNeighbor = v->getNeighborNum();
        // Voronoi cell area
        double exDiff = 0.0;
        exDiff += 2 * dotProduct(
            v->force,
            array2Vector<double, 3>(v->getMVoronoiCell()->getDArea())
        );
        for(size_t vIdx = 0; vIdx < numNeighbor; ++vIdx) {
            exDiff += 2 * dotProduct(
                v->getNeighborVertices()[vIdx]->force,
                array2Vector<double, 3>(v->getMVoronoiCell()->getDNeighborArea()[vIdx])
            );
        }
        EXPECT_NEAR(vCellArea1[idx] - vCellArea2[idx], exDiff, abs(exDiff / 1000));
        // Voronoi cell curvature
        exDiff = 0.0;
        exDiff += 2 * dotProduct(
            v->force,
            array2Vector<double, 3>(v->getMVoronoiCell()->getDCurv())
        );
        for(size_t vIdx = 0; vIdx < numNeighbor; ++vIdx) {
            exDiff += 2 * dotProduct(
                v->getNeighborVertices()[vIdx]->force,
                array2Vector<double, 3>(v->getMVoronoiCell()->getDNeighborCurv()[vIdx])
            );
        }
        EXPECT_NEAR(vCellCurv1[idx] - vCellCurv2[idx], exDiff, abs(exDiff / 1000));
    }

}


#  endif //DO_THIS_GEOMETRY_MEMBRANE_TEST
#endif //TESTING

