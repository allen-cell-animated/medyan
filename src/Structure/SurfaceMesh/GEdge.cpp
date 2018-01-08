#include "Edge.h"
#include "GEdge.h"

#include "Triangle.h"
#include "GTriangle.h"

#include "MathFunctions.h"

using namespace mathfunc;

void GEdge::calcLength() {
    auto& v0 = _pEdge->getVertices()[0];
    auto& v1 = _pEdge->getVertices()[1];

    _currentLength = twoPointDistance(v0->coordinate, v1->coordinate);
    for(int coordIdx = 0; coordIdx < 3; ++coordIdx) {
        _dCurrentLength[0][coordIdx] = (v0->coordinate[coordIdx] - v1->coordinate[coordIdx]) / _currentLength;
        _dCurrentLength[1][coordIdx] = (v1->coordinate[coordIdx] - v0->coordinate[coordIdx]) / _currentLength;
    }
}

void GEdge::calcStretchedLength(double d) {
    auto& v0 = _pEdge->getVertices()[0];
    auto& v1 = _pEdge->getVertices()[1];
    _stretchedLength = twoPointDistanceStretched(v0->coordinate, v0->force,
                                                 v1->coordinate, v1->force,
                                                 d);
}

void GEdge::calcPseudoUnitNormal() {
    /*
    This function depends on the result of triangle unit normal
    */
    GTriangle* gt0 = _pEdge->getTriangles()[0]->getGTriangle();
    GTriangle* gt1 = _pEdge->getTriangles()[1]->getGTriangle();
    _pseudoUnitNormal = vectorSum(gt0->getUnitNormal(), gt1->getUnitNormal());
    normalize(_pseudoUnitNormal);
}

void GEdge::calcStretchedPseudoUnitNormal(double d) {
    /*
    This function depends on the result of stretched triangle unit normal
    */
    GTriangle* gt0 = _pEdge->getTriangles()[0]->getGTriangle();
    GTriangle* gt1 = _pEdge->getTriangles()[1]->getGTriangle();
    _stretchedPseudoUnitNormal = vectorSum(gt0->getStretchedUnitNormal(), gt1->getStretchedUnitNormal());
    normalize(_stretchedPseudoUnitNormal);
}
