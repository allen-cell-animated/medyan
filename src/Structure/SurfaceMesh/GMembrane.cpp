#include "GMembrane.h"

#include "Membrane.h"
#include "Triangle.h"
#include "MTriangle.h"
#include "Vertex.h"

void GMembrane::calcVolume() {
    auto& vs = _pMembrane->getVertexVector();
    size_t numVertices = vs.size();

    _dVolume.resize(numVertices);
    _Volume = 0;
    std::fill(_dVolume.begin(), _dVolume.end(), 0.0);

    for(Triangle* t: _pMembrane->getTriangleVector()) {
        // TODO:
    }
}
