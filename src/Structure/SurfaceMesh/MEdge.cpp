#include "MEdge.h"

#include "MathFunctions.h"

void MEdge::calcLength() {
    auto& b0 = _he[0]->_b[0];
    auto& b1 = _he[0]->_b[1];

    _currentLength = twoPointDistance(b0->coordinate, b1->coordinate);
    for(int coordIdx = 0; coordIdx < 2; ++coordIdx) {
        _dCurrentLength[0][coordIdx] = (b0->coordinate[coordIdx] - b1->coordinate[coordIdx]) / _currentLength;
        _dCurrentLength[1][coordIdx] = (b1->coordinate[coordIdx] - b0->coordinate[coordIdx]) / _currentLength;
    }
}