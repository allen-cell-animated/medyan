#ifndef MEDYAN_GVoronoiCell_h
#define MEDYAN_GVoronoiCell_h

#include <array>
#include <vector>

class Vertex; // Forward declaration

class GVoronoiCell {

private:
    Vertex* _pVertex; // Parent vertex

    // Note: vectors must initially be resized to the number of neighbors around the vertex.
    double _currentArea; // Current area
    std::array<double, 3> _dCurrentArea; // Derivative of area on the central vertex
    std::vector<std::array<double, 3>> _dNeighborCurrentArea; // Derivative of the area on the neighboring vertices
    double _stretchedArea; // Temporarily store the stretched area

    double _currentCurv; // Current mean curvature
    std::array<double, 3> _dCurrentCurv;
    std::vector<std::array<double, 3>> _dNeighborCurrentCurv; // Derivative of the mean curv on the neighboring vertices
    double _stretchedCurv; // Temporarily store the stretched mean curvature

    std::array<double, 3> _pseudoUnitNormal; // Pseudo unit normal around the vertex
    std::array<double, 3> _stretchedPseudoUnitNormal; // Pseudo unit normal under stretch

public:
    GVoronoiCell(size_t numNeighbors);

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
