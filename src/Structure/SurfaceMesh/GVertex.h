#ifndef MEDYAN_GVertex_h
#define MEDYAN_GVertex_h

#include <array>
#include <vector>

#include "MathFunctions.h"

class Vertex; // Forward declaration

struct GVertex {

    double area; // Current area
    mathfunc::Vec3 dArea; // Derivative of area on the central vertex, derivative on neighbors are stored in half edges
    double sArea; // Temporarily store the stretched area

    double curv; // Current mean curvature
    mathfunc::Vec3 dCurv;
    double sCurv; // Temporarily store the stretched mean curvature

    mathfunc::Vec3 dVolume; // Derivative of volume on this vertex

    mathfunc::Vec3 pseudoUnitNormal; // Pseudo unit normal around the vertex
    mathfunc::Vec3 sPseudoUnitNormal; // Pseudo unit normal under stretch


    /// Set parent 
    void setVertex(Vertex* v) { _pVertex = v; }
    Vertex* getVertex() { return _pVertex; }

    double getArea() { return _currentArea; }
    std::array<double, 3>& getDArea() { return _dCurrentArea; }
    std::vector<std::array<double, 3>>& getDNeighborArea() { return _dNeighborCurrentArea; }
    void calcArea();
    double getStretchedArea() { return _stretchedArea; }
    void calcStretchedArea(double d); // Calculates the stretched area, and store the result in _stretchedArea.
                                      // Does not calculate the derivatives.

    double getCurv() { return _currentCurv; }
    std::array<double, 3>& getDCurv() { return _dCurrentCurv; }
    std::vector<std::array<double, 3>>& getDNeighborCurv() { return _dNeighborCurrentCurv; }
    void calcCurv();
    double getStretchedCurv() { return _stretchedCurv; }
    void calcStretchedCurv(double d); // Calculates the stretched mean curvature, and store the result in _stretchedCurv.
                                      // Does not calculate the derivatives.
    
    std::array<double, 3>& getPseudoUnitNormal() { return _pseudoUnitNormal; }
    void calcPseudoUnitNormal(); // Calculates the pseudo unit normal w/o derivatives
    std::array<double, 3>& getStretchedPseudoUnitNormal() { return _stretchedPseudoUnitNormal; }
    void calcStretchedPseudoUnitNormal(double d); // Calculates the pseudo unit normal w/o derivatives

};


#endif
