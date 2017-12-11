#include <algorithm>
#include <functional>

#include "Vertex.h"
#include "MVoronoiCell.h"

#include "MathFunctions.h"

using namespace mathfunc;

void MVoronoiCell::calcArea() {
    /**************************************************************************
    This calculation depends on the result of
        - the angle calculation of triangles
        - the length calculation of edges

    Currently only Voronoi area is considered, but we can also change it to
    "Mixed area" (Meyer 2003). Voronoi area can become 0 or even negative when
    obtuse triangles are encountered.

    The area around vertex i is:
    A = 1/8 * \sum_{v_j as v_i's neighbor} (cot \alpha_ij + cot \beta_ij ) * |r_i - r_j|^2
    (See Meyer 2003). As a result, the area is dependent on the center vertex,
    as well as all the neighboring vertices tethered together.
    **************************************************************************/
    
    size_t n = _pVertex -> getNeighborNum();
    _currentArea = 0;
    for(size_t nIdx = 0; nIdx < n; ++nIdx){
        for(int coordIdx = 0; coordIdx < 3; ++coordIdx){
            _dNeighborCurrentArea[nIdx][coordIdx] = 0;
        }
    }

    for(size_t nIdx = 0; nIdx < n; ++nIdx){
        // Think about the triangle (center, nIdx, nIdx+1)
        auto mTriangleR = _pVertex->getNeighborTriangles()[nIdx]->getMTriangle();
        auto mTriangleL = _pVertex->getNeighborTriangles()[(nIdx + n - 1) % n]->getMTriangle();
        int triRIdx0 = (3 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx1 = (4 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx2 = (5 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triLIdx0 = (3 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx1 = (4 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx2 = (5 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;

        // Think about the edge (center, nIdx)
        auto mEdge = _pVertex->getNeighborEdges()[nIdx]->getMEdge();
        int edgeIdx0 = (2 - _pVertex->getEdgeHead()[nIdx]) % 2;
        int edgeIdx1 = (3 - _pVertex->getEdgeHead()[nIdx]) % 2;

        double dist2 = mEdge->getLength() * mEdge->getLength();
        double sumCotTheta = mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2];
        _currentArea += sumCotTheta * dist2;
        for(int coordIdx = 0; coordIdx < 3; ++coordIdx){
            _dCurrentArea[coordIdx] += ((mTriangleL->getDCotTheta()[triLIdx1][triLIdx0][coordIdx] + mTriangleR->getDCotTheta()[triRIdx2][triRIdx0][coordIdx]) * dist2
                + sumCotTheta * 2 * mEdge->getLength() * mEdge->getDLength()[edgeIdx0][coordIdx]) / 8;
            _dNeighborCurrentArea[nIdx][coordIdx] += ((mTriangleL->getDCotTheta()[triLIdx1][triLIdx2][coordIdx] + mTriangleR->getDCotTheta()[triRIdx2][triRIdx1][coordIdx]) * dist2
                + sumCotTheta * 2 * mEdge->getLength() * mEdge->getDLength()[edgeIdx1][coordIdx]) / 8;
            _dNeighborCurrentArea[(nIdx + n - 1) % n][coordIdx] += mTriangleL->getDCotTheta()[triLIdx1][triLIdx1][coordIdx] * dist2 / 8;
            _dNeighborCurrentArea[(nIdx + 1) % n][coordIdx] += mTriangleR->getDCotTheta()[triRIdx2][triRIdx2][coordIdx] * dist2 / 8;
        }
        
    }
}

void MVoronoiCell::calcStretchedArea(double d) {
    /**************************************************************************
    This calculation depends on the result of
        - the stretched angle calculation of triangles
        - the stretched length calculation of edges

    The formulae used here must match those used in calcArea()

    Coincidentally, variable d is not used here, because the information could
    be obtained from edges and triangles. However, for consistency and possible
    modification in the future, this variable is still listed here as a
    function parameter.
    **************************************************************************/
    
    size_t n = _pVertex -> getNeighborNum();
    _stretchedArea = 0;

    for(size_t nIdx = 0; nIdx < n; ++nIdx){
        // Think about the triangle (center, nIdx, nIdx+1)
        // The commented indices are not used, but serve as a reference of what their relations are.
        auto mTriangleR = _pVertex->getNeighborTriangles()[nIdx]->getMTriangle();
        auto mTriangleL = _pVertex->getNeighborTriangles()[(nIdx + n - 1) % n]->getMTriangle();
        //int triRIdx0 = (3 - _pVertex->getTriangleHead()[nIdx]) % 3;
        //int triRIdx1 = (4 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx2 = (5 - _pVertex->getTriangleHead()[nIdx]) % 3;
        //int triLIdx0 = (3 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx1 = (4 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        //int triLIdx2 = (5 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;

        // Think about the edge (center, nIdx)
        auto mEdge = _pVertex->getNeighborEdges()[nIdx]->getMEdge();

        double dist2Stretched = mEdge->getStretchedLength() * mEdge->getStretchedLength();
        double sumCotThetaStretched = mTriangleL->getStretchedCotTheta()[triLIdx1] + mTriangleR->getStretchedCotTheta()[triRIdx2];
        _stretchedArea += sumCotThetaStretched * dist2Stretched;        
    }
}

void MVoronoiCell::calcCurv() {
    /**************************************************************************
    This calculation depends on the result of
        - the area calculation of the current Voronoi cell
        - the angle calculation of triangles
        //- the length calculation of edges

    The discretized Laplace-Beltrami operator operated on the local vertex is
    K = 2 * Curv * n, where n is the normal vector, and
    K = 1/(2A) * \sum_{v_j as v_i's neighbor} (cot \alpha_ij + cot \beta_ij ) * (r_i - r_j)
    (See Meyer 2003). As a result, the mean curvature is dependent on the
    center vertex, as well as all the neighboring vertices tethered together.
    **************************************************************************/

    size_t n = _pVertex -> getNeighborNum();

    std::array<double, 3> k = {}; // Result of Laplace-Beltrami operator
    std::array<std::array<double, 3>, 3> dK = {};
    std::vector<std::array<std::array<double, 3>, 3>> dNeighborK(n, {{}}); // Use {{}} instead of {} here to ensure it's not an allocator type
                                                                           // For forward compatibility starting from C++14

    for(size_t nIdx = 0; nIdx < n; ++nIdx){
        // Think about the triangle (center, nIdx, nIdx+1)
        auto mTriangleR = _pVertex->getNeighborTriangles()[nIdx]->getMTriangle();
        auto mTriangleL = _pVertex->getNeighborTriangles()[(nIdx + n - 1) % n]->getMTriangle();
        int triRIdx0 = (3 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx1 = (4 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx2 = (5 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triLIdx0 = (3 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx1 = (4 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx2 = (5 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;

        std::array<double, 3> diff;
        double sumCotTheta = mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2];
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            diff[coordIdx] = _pVertex->coordinate[coordIdx] - _pVertex->getNeighborVertices()[nIdx]->coordinate[coordIdx];
            k[coordIdx] += sumCotTheta * diff[coordIdx];
        }

        // Now, dDiff is Eye3, and dNDiff is -Eye3
        std::array<std::array<double, 3>, 3> tensorTmp0 = tensorProduct<3>(vectorSum<3>(mTriangleL->getDCotTheta()[triLIdx1][triLIdx0], mTriangleR->getDCotTheta()[triRIdx2][triRIdx0]), diff);
        std::array<std::array<double, 3>, 3> tensorTmp1 = tensorProduct<3>(vectorSum<3>(mTriangleL->getDCotTheta()[triLIdx1][triLIdx2], mTriangleR->getDCotTheta()[triRIdx2][triRIdx1]), diff);
        std::array<std::array<double, 3>, 3> tensorTmpL = tensorProduct<3>(mTriangleL->getDCotTheta()[triLIdx1][triLIdx1], diff);
        std::array<std::array<double, 3>, 3> tensorTmpR = tensorProduct<3>(mTriangleR->getDCotTheta()[triRIdx2][triRIdx2], diff);

        matrixIncrease<3>(dK, matrixSum<3>(tensorTmp0, matrixMultiply<3>(Eye3, sumCotTheta)));
        matrixIncrease<3>(dNeighborK[nIdx], matrixDifference<3>(tensorTmp1, matrixMultiply<3>(Eye3, sumCotTheta)));
        matrixIncrease<3>(dNeighborK[(nIdx + n - 1) % n], tensorTmpL);
        matrixIncrease<3>(dNeighborK[(nIdx + 1) % n], tensorTmpR);

    }

    // Convert K to K/2A
    dK = matrixMultiply<3>(matrixDifference<3>(matrixMultiply<3>(dK, _currentArea), tensorProduct<3>(_dCurrentArea, k)), 0.5 / _currentArea / _currentArea);
    for(size_t nIdx = 0; nIdx < n; ++nIdx) {
        dNeighborK[nIdx] = matrixMultiply<3>(matrixDifference<3>(matrixMultiply<3>(dNeighborK[nIdx], _currentArea), tensorProduct<3>(_dNeighborCurrentArea[nIdx], k)), 0.5 / _currentArea / _currentArea);
    }
    vectorExpand<3>(k, 0.5 / _currentArea);

    // Calculate mean curvature
    _currentCurv = magnitude<3>(k) / 2;
    _dCurrentCurv = vectorMultiply<3>(matrixProduct<3>(dK, k), 1.0 / 4 / _currentCurv);
    for(size_t nIdx = 0; nIdx < n; ++nIdx) {
        _dNeighborCurrentCurv[nIdx] = vectorMultiply<3>(matrixProduct<3>(dNeighborK[nIdx], k), 1.0 / 4 / _currentCurv);
    }
}

void MVoronoiCell::calcStretchedCurv(double d) {
    /**************************************************************************
    This calculation depends on the result of
        - the stretched area calculation of the current Voronoi cell
        - the stretched angle calculation of triangles

    The formulae used here must match those used in calcCurv()
    **************************************************************************/

    size_t n = _pVertex -> getNeighborNum();

    std::array<double, 3> kStretched = {}; // Result of Laplace-Beltrami operator

    for(size_t nIdx = 0; nIdx < n; ++nIdx){
        // Think about the triangle (center, nIdx, nIdx+1)
        // The commented indices are not used, but serve as a reference of what their relations are.
        auto mTriangleR = _pVertex->getNeighborTriangles()[nIdx]->getMTriangle();
        auto mTriangleL = _pVertex->getNeighborTriangles()[(nIdx + n - 1) % n]->getMTriangle();
        //int triRIdx0 = (3 - _pVertex->getTriangleHead()[nIdx]) % 3;
        //int triRIdx1 = (4 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx2 = (5 - _pVertex->getTriangleHead()[nIdx]) % 3;
        //int triLIdx0 = (3 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx1 = (4 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        //int triLIdx2 = (5 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;

        std::array<double, 3> diffStretched;
        double sumCotThetaStretched = mTriangleL->getStretchedCotTheta()[triLIdx1] + mTriangleR->getStretchedCotTheta()[triRIdx2];
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            diffStretched[coordIdx] = _pVertex->coordinate[coordIdx] + d * _pVertex->force[coordIdx]
                                    - (_pVertex->getNeighborVertices()[nIdx]->coordinate[coordIdx] + d * _pVertex->getNeighborVertices()[nIdx]->force[coordIdx]);
            kStretched[coordIdx] += sumCotThetaStretched * diffStretched[coordIdx];
        }
    }

    // Convert K to K/2A
    vectorExpand<3>(kStretched, 0.5 / _stretchedArea);

    // Calculate mean curvature
    _stretchedCurv = magnitude<3>(kStretched) / 2;
}
