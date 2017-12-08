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
    for(int nIdx = 0; nIdx < n; ++nIdx){
        for(int coordIdx = 0; coordIdx < 3; ++coordIdx){
            _dNeighborCurrentArea[nIdx][coordIdx] = 0;
        }
    }

    for(int nIdx = 0; nIdx < n; ++nIdx){
        // Think about the triangle (center, nIdx, nIdx+1)
        auto& mTriangleR = _pVertex->getNeighborTriangles()[nIdx]->getMTriangle();
        auto& mTriangleL = _pVertex->getNeighborTriangles()[(nIdx + n - 1) % n]->getMTriangle();
        int triRIdx0 = (3 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx1 = (4 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx2 = (5 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triLIdx0 = (3 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx1 = (4 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx2 = (5 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        // Think about the edge (center, nIdx)
        auto& mEdge = _pVertex->getNeighborEdges()[nIdx]->getMEdge();
        int edgeIdx0 = (2 - _pVertex->getEdgeHead()[nIdx]) % 2;
        int edgeIdx1 = (3 - _pVertex->getEdgeHead()[nIdx]) % 2;

        double dist2 = mEdge->getLength() * mEdge->getLength();
        double sumCotTheta = mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2];
        _currentArea += sumCotTheta * dist2;
        for(int coordIdx = 0; coordIdx < 3; ++coordIdx){
            _dCurrentArea += ((mTriangleL->getDCotTheta()[triLIdx1][triLIdx0][coordIdx] + mTriangleR->getDCotTheta()[triRIdx2][triRIdx0][coordIdx]) * dist2
                + sumCotTheta * 2 * mEdge->getLength() * mEdge->getDLength()[edgeIdx0][coordIdx]) / 8;
            _dNeighborCurrentArea[nIdx][coordIdx] += ((mTriangleL->getDCotTheta()[triLIdx1][triLIdx2][coordIdx] + mTriangleR->getDCotTheta()[triRIdx2][triRIdx1][coordIdx]) * dist2
                + sumCotTheta * 2 * mEdge->getLength() * mEdge->getDLength()[edgeIdx1][coordIdx]) / 8;
            _dNeighborCurrentArea[(nIdx + n - 1) % n][coordIdx] += mTriangleL->getDCotTheta()[triLIdx1][triLIdx1] * dist2 / 8;
            _dNeighborCurrentArea[(nIdx + 1) % n][coordIdx] += mTriangleR->getDCotTheta()[triRIdx2][triRIdx2] * dist2 / 8;
        }
        
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

    // Use vector instead of array here for tensor product templates.
    std::vector<double> k(3, 0.0); // Result of Laplace-Beltrami operator
    std::vector<std::vector<double>> dK(3, std::vector<double>(3, 0.0));
    std::vector<std::vector<std::vector<double>>> dNeighborK(n, std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0)));

    for(size_t nIdx = 0; nIdx < n; ++nIdx){
        // Think about the triangle (center, nIdx, nIdx+1)
        auto& mTriangleR = _pVertex->getNeighborTriangles()[nIdx]->getMTriangle();
        auto& mTriangleL = _pVertex->getNeighborTriangles()[(nIdx + n - 1) % n]->getMTriangle();
        int triRIdx0 = (3 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx1 = (4 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triRIdx2 = (5 - _pVertex->getTriangleHead()[nIdx]) % 3;
        int triLIdx0 = (3 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx1 = (4 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        int triLIdx2 = (5 - _pVertex->getTriangleHead()[(nIdx + n - 1) % n]) % 3;
        // Think about the edge (center, nIdx)
        auto& mEdge = _pVertex->getNeighborEdges()[nIdx]->getMEdge();
        int edgeIdx0 = (2 - _pVertex->getEdgeHead()[nIdx]) % 2;
        int edgeIdx1 = (3 - _pVertex->getEdgeHead()[nIdx]) % 2;

        std::vector<double> diff(3);
        double sumCotTheta = mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2];
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            diff[coordIdx] = _pVertex->coordinate[coordIdx] - _pVertex->getNeighborVertices()[nIdx]->coordinate[coordIdx];
            k[coordIdx] += sumCotTheta * diff[coordIdx];
        }

        // Now, dDiff is Eye3, and dNDiff is -Eye3
        std::vector<std::vector<double>> tensorTmp0 = tensorProduct(vectorSum(mTriangleL->getDCotTheta[triLIdx1][triLIdx0], mTriangleR->getDCotTheta[triRIdx2][triRIdx0]), diff);
        std::vector<std::vector<double>> tensorTmp1 = tensorProduct(vectorSum(mTriangleL->getDCotTheta[triLIdx1][triLIdx2], mTriangleR->getDCotTheta[triRIdx2][triRIdx1]), diff);
        std::vector<std::vector<double>> tensorTmpL = tensorProduct(mTriangleL->getDCotTheta[triLIdx1][triLIdx1], diff);
        std::vector<std::vector<double>> tensorTmpR = tensorProduct(mTriangleR->getDCotTheta[triRIdx2][triRIdx2], diff);

        matrixIncrease(dK, matrixSum(tensorTmp0, matrixMultiply(Eye3, sumCotTheta)));
        matrixIncrease(dNeighborK[nIdx], matrixDifference(tensorTmp1, matrixMultiply(Eye3, sumCotTheta)));
        matrixIncrease(dNeighborK[(nIdx + n - 1) % n], tensorTmpL);
        matrixIncrease(dNeighborK[(nIdx + 1) % n], tensorTmpR);

    }

    // Convert K to K/2A
    dK = matrixMultiply(matrixDifference(matrixMultiply(dK, _currentArea), tensorProduct(array2Vector<double, 3>(_dCurrentArea), k)), 0.5 / _currentArea / _currentArea);
    for(size_t nIdx = 0; nIdx < n; ++nIdx) {
        dNeighborK[nIdx] = matrixMultiply(matrixDifference(matrixMultiply(dNeighborK[nIdx], _currentArea), tensorProduct(array2Vector<double, 3>(_dNeighborCurrentArea[nIdx]), k)), 0.5 / _currentArea / _currentArea);
    }
    vectorExpand(k, 0.5 / _currentArea);

    // Calculate mean curvature
    _currentCurv = magnitude(k) / 2;
    _dCurrentCurv = vectorMultiply(matrixProduct(dK, k), 1.0 / 4 / _currentCurv);
    for(size_t nIdx = 0; nIdx < n; ++nIdx) {
        _dNeighborCurrentCurv[nIdx] = vectorMultiply(matrixProduct(dNeighborK[nIdx], k), 1.0 / 4 / _currentCurv);
    }
}
