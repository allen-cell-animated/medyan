
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

TEST_F(MembraneGeometryTest, CheckGeometry) {

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
    double exVCellCurv = 1;
    for(Vertex* it: m->getVertexVector()) {
        EXPECT_DOUBLE_EQ(it->getMVoronoiCell()->getArea(), exVCellArea);
        cout << it->getMVoronoiCell()->getCurv();
    }
}


#  endif //DO_THIS_GEOMETRY_MEMBRANE_TEST
#endif //TESTING

