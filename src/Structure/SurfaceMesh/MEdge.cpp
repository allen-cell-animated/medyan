#include "MEdge.h"

#include "MathFunctions.h"

void MEdge::calcLength() {
    auto& b0 = _pEdge->_v[0]->_b;
    auto& b1 = _pEdge->_v[1]->_b;

    _currentLength = twoPointDistance(b0->coordinate, b1->coordinate);
    for(int coordIdx = 0; coordIdx < 3; ++coordIdx) {
        _dCurrentLength[0][coordIdx] = (b0->coordinate[coordIdx] - b1->coordinate[coordIdx]) / _currentLength;
        _dCurrentLength[1][coordIdx] = (b1->coordinate[coordIdx] - b0->coordinate[coordIdx]) / _currentLength;
    }
}