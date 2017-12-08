#include "Edge.h"
#include "MEdge.h"

#include "MathFunctions.h"

using namespace mathfunc;

void MEdge::calcLength() {
    auto& v0 = _pEdge->getVertices()[0];
    auto& v1 = _pEdge->getVertices()[1];

    _currentLength = twoPointDistance(v0->coordinate, v1->coordinate);
    for(int coordIdx = 0; coordIdx < 3; ++coordIdx) {
        _dCurrentLength[0][coordIdx] = (v0->coordinate[coordIdx] - v1->coordinate[coordIdx]) / _currentLength;
        _dCurrentLength[1][coordIdx] = (v1->coordinate[coordIdx] - v0->coordinate[coordIdx]) / _currentLength;
    }
}