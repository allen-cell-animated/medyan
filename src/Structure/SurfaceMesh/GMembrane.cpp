#include "GMembrane.h"

#include "MathFunctions.h"
using namespace mathfunc;

#include "Membrane.h"
#include "Triangle.h"
#include "MTriangle.h"
#include "Vertex.h"

void GMembrane::calcVolume() {
    auto& vs = _pMembrane->getVertexVector();
    size_t numVertices = vs.size();

    _dVolume.resize(numVertices);
    _volume = 0.0;
	std::fill(_dVolume.begin(), _dVolume.end(), array<double, 3>{});

    for(Triangle* t: _pMembrane->getTriangleVector()) {
        auto v0 = vector2Array<double, 3>(t->getVertices()[0]->coordinate);
        auto v1 = vector2Array<double, 3>(t->getVertices()[1]->coordinate);
        auto v2 = vector2Array<double, 3>(t->getVertices()[2]->coordinate);

        static constexpr double oneOverSix = 1.0 / 6.0;

        auto r0 = vectorMultiply(crossProduct(v1, v2), oneOverSix);
        auto r1 = vectorMultiply(crossProduct(v2, v0), oneOverSix);
        auto r2 = vectorMultiply(crossProduct(v0, v1), oneOverSix);

        _volume += dotProduct(v0, r0);

        vectorIncrease(_dVolume[t->getVertices()[0]->getMembraneVertexIdx()], r0);
        vectorIncrease(_dVolume[t->getVertices()[1]->getMembraneVertexIdx()], r1);
        vectorIncrease(_dVolume[t->getVertices()[2]->getMembraneVertexIdx()], r2);
    }
}

void GMembrane::calcStretchedVolume(double d) {
    _stretchedVolume = 0.0;

    for(Triangle* t: _pMembrane->getTriangleVector()) {
        auto v0 = vectorSum(
            vector2Array<double, 3>(t->getVertices()[0]->coordinate),
            vectorMultiply(vector2Array<double, 3>(t->getVertices()[0]->force), d)
        );
        auto v1 = vectorSum(
            vector2Array<double, 3>(t->getVertices()[1]->coordinate),
            vectorMultiply(vector2Array<double, 3>(t->getVertices()[1]->force), d)
        );
        auto v2 = vectorSum(
            vector2Array<double, 3>(t->getVertices()[2]->coordinate),
            vectorMultiply(vector2Array<double, 3>(t->getVertices()[2]->force), d)
        );

        static constexpr double oneOverSix = 1.0 / 6.0;

        _stretchedVolume += dotProduct(v0, crossProduct(v1, v2)) * oneOverSix;
    }
}
