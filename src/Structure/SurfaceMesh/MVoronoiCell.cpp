#include "MVoronoiCell.h"

#include "MathFunctions.h"

void MVoronoiCell::calcArea() {
    /**************************************************************************
    This calculation depends on the result of the angle calculation of triangles.

    Currently only Voronoi area is considered, but we can also change it to
    "Mixed area" (Meyer 2003). Voronoi area can become 0 or even negative when
    obtuse triangles are encountered.

    The area around vertex i is:
    A = 1/8 * \sum_{v_j as v_i's neighbor} (cot \alpha_ij + cot \beta_ij ) * |r_i - r_j|^2
    (See Meyer 2003). As a result, the area is dependent on the center vertex,
    as well as all the neighboring vertices tethered together.
    **************************************************************************/
    
    _currentArea = 0;

    // TODO: implement this
    
    auto& v0 = _pTriangle->getVertices()[0]->getBead()->coordinate;
    auto& v1 = _pTriangle->getVertices()[1]->getBead()->coordinate;
    auto& v2 = _pTriangle->getVertices()[2]->getBead()->coordinate;

    double l0 = _pTriangle->getEdges()[0]->getMEdge()->getLength(); // length v0 - v1
    double l1 = _pTriangle->getEdges()[1]->getMEdge()->getLength(); // length v1 - v2
    double l2 = _pTriangle->getEdges()[2]->getMEdge()->getLength(); // length v2 - v0

    double dot12 = scalarProduct(v0, v1, v0, v2);

    _currentArea = 0.5 * magnitude(vectorProduct(v0, v1, v0, v2));

    // Calculate gradients
    for(int coordIdx = 0; coordIdx < 3; ++coordIdx) {
        double tmpR01 = v1[coordIdx] - v0[coordIdx];
        double tmpR02 = v2[coordIdx] - v0[coordIdx];

        _dCurrentArea[0][coordIdx] = (-l0*l0*tmpR02 - l2*l2*tmpR01 + dot12*(tmpR01 + tmpR02)) / _currentArea / 4;
        _dCurrentArea[1][coordIdx] = (l2*l2*tmpR01 - dot12*tmpR02) / _currentArea / 4;
        _dCurrentArea[2][coordIdx] = (l0*l0*tmpR02 - dot12*tmpR01) / _currentArea / 4;
    }

}
