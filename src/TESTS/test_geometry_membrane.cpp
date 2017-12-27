
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

#    include "SubSystem.h"

#    include "Vertex.h"
#    include "Edge.h"
#    include "Triangle.h"
#    include "MVoronoiCell.h"
#    include "MEdge.h"
#    include "MTriangle.h"
#    include "Membrane.h"

using vertexData = tuple<array<double, 3>, vector<size_t>>;
using membraneData = vector<vertexData>;

membraneData membraneDataOctahedron(double radius) {
    return membraneData{
        vertexData({0, 0, radius}, {1, 2, 3, 4}),
        vertexData({radius, 0, 0}, {0, 4, 5, 2}),
        vertexData({0, radius, 0}, {0, 1, 5, 3}),
        vertexData({-radius, 0, 0}, {0, 2, 5, 4}),
        vertexData({0, -radius, 0}, {0, 3, 5, 1}),
        vertexData({0, 0, -radius}, {4, 3, 2, 1})
    };
}

TEST(TopologyTest, Main) {
    SubSystem s;
    membraneData memData = membraneDataOctahedron(100);
    Membrane m(&s, 0, memData);

    // Check that the vertices, edges and triangles are correctly registered.
    EXPECT_EQ(m.getVertexVector().size(), 6);
    EXPECT_EQ(m.getEdgeVector().size(), 12);
    EXPECT_EQ(m.getTriangleVector().size(), 8);

    
}


#  endif //DO_THIS_GEOMETRY_MEMBRANE_TEST
#endif //TESTING

