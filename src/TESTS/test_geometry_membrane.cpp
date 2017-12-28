
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

#    include "common.h"

#    include "GController.h"
#    include "SubSystem.h"

#    include "Vertex.h"
#    include "Edge.h"
#    include "Triangle.h"
#    include "MVoronoiCell.h"
#    include "MEdge.h"
#    include "MTriangle.h"
#    include "Membrane.h"

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
public:
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
        EXPECT_DOUBLE_EQ(it->getMEdge()->getLength(), exEdgeLen);
    // Check triangle area and angle
    double exTriangleArea = radius * radius * sqrt(3) / 2;
    double exTriangleAngle = M_PI / 3;
    for(Triangle* it: m->getTriangleVector()) {
        EXPECT_DOUBLE_EQ(it->getMTriangle()->getArea(), exTriangleArea);
        for(double eachTheta: it->getMTriangle()->getTheta())
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

    m->updateGeometry(false, 1.0);

    // Check edge length
    double exStretchedEdgeLen = radius;
    for(Edge* it: v0->getNeighborEdges())
        EXPECT_DOUBLE_EQ(it->getMEdge()->getStretchedLength(), exStretchedEdgeLen);
    // Check triangle area and angle
    double exStretchedTriangleArea = radius * radius / 2;
    double exStretchedTriangleAngleIn = M_PI / 2;
    double exStretchedTriangleAngleOut = M_PI / 4;
    for(size_t nIdx = 0; nIdx < numNeighborV0; ++nIdx) {
        EXPECT_DOUBLE_EQ(v0->getNeighborTriangles()[nIdx]->getMTriangle()->getStretchedArea(), exStretchedTriangleArea);
        for(size_t angleIdx = 0; angleIdx < 3; ++angleIdx) {
            if((3 - v0->getTriangleHead()[nIdx]) % 3 == angleIdx)
                EXPECT_DOUBLE_EQ(
                    v0->getNeighborTriangles()[nIdx]->getMTriangle()->getStretchedTheta()[angleIdx],
                    exStretchedTriangleAngleIn
                );
            else
                EXPECT_DOUBLE_EQ(
                    v0->getNeighborTriangles()[nIdx]->getMTriangle()->getStretchedTheta()[angleIdx],
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
}


#  endif //DO_THIS_GEOMETRY_MEMBRANE_TEST
#endif //TESTING

