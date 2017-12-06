#include "MVoronoiCell.h"

#include "MathFunctions.h"

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
    
    n = _pVertex -> getNeighborNum();
    _currentArea = 0;
    for(int nIdx = 0; nIdx < n; ++nIdx){
        for(int coordIdx = 0; coordIdx < 3; ++coordIdx){
            _dTetheredCurrentArea[nIdx][coordIdx] = 0;
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
        _currentArea += (mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2]) * dist2;
        for(int coordIdx = 0; coordIdx < 3; ++coordIdx){
            _dCurrentArea += ((mTriangleL->getDCotTheta()[triLIdx1][triLIdx0][coordIdx] + mTriangleR->getDCotTheta()[triRIdx2][triRIdx0][coordIdx]) * dist2
                + (mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2]) * 2 * mEdge->getLength() * mEdge->getDLength()[edgeIdx0][coordIdx]) / 8;
            _dTetheredCurrentArea[nIdx][coordIdx] += ((mTriangleL->getDCotTheta()[triLIdx1][triLIdx2][coordIdx] + mTriangleR->getDCotTheta()[triRIdx2][triRIdx1][coordIdx]) * dist2
                + (mTriangleL->getCotTheta()[triLIdx1] + mTriangleR->getCotTheta()[triRIdx2]) * 2 * mEdge->getLength() * mEdge->getDLength()[edgeIdx1][coordIdx]) / 8;
            _dTetheredCurrentArea[(nIdx + n - 1) % n][coordIdx] += mTriangleL->getDCotTheta()[triLIdx1][triLIdx1] * dist2 / 8;
            _dTetheredCurrentArea[(nIdx + 1) % n][coordIdx] += mTriangleR->getDCotTheta()[triRIdx2][triRIdx2] * dist2 / 8;
        }
        
    }
}

void MVoronoiCell::calcCurv() {

    // TODO: Implement this
}