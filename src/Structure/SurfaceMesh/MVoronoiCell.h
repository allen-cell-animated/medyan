#ifndef MEDYAN_MVoronoiCell_h
#define MEDYAN_MVoronoiCell_h

#include <array>
#include <vector>

class Vertex; // Forward declaration

class MVoronoiCell {

private:
    Vertex* _pVertex; // Parent vertex

    double _eqArea; // Length of unstretched area, determined at mesh generation
    double _kElastic; // Local elastic modulus
    double _kBending; // Local bending modulus
    double _eqCurv; // Local spontaneous curvature

    // Note: vectors must 
    double _currentArea; // Current area
    std::array<double, 3> _dCurrentArea; // Derivative of area on the central vertex
    std::vector<std::array<double, 3>> _dNeighborCurrentArea; // Derivative of the area on the neighboring vertices
    double _stretchedArea; // Temporarily store the stretched area

    double _currentCurv; // Current mean curvature
    std::array<double, 3> _dCurrentCurv;
    std::vector<std::array<double, 3>> _dNeighborCurrentCurv; // Derivative of the mean curv on the neighboring vertices
    double _stretchedCurv; // Temporarily store the stretched mean curvature

public:
    MVoronoiCell(double eqArea) {
        setEqArea(eqArea);
    }

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setElasticModulus(double kElastic) { _kElastic = kElastic; }
    double getElasticModulus() { return _kElastic; }

    void setBendingModulus(double kBending) { _kBending = kBending; }
    double getBendingModulus() { return _kBending; }

    void setEqCurv(double eqCurv) { _eqCurv = eqCurv; }
    double getEqCurv() { return _eqCurv; }

    double getArea() { return _currentArea; }
    void calcArea();
    double getStretchedArea() { return _stretchedArea; }
    void calcStretchedArea(double d); // Calculates the stretched area, and store the result in _stretchedArea.
                                      // Does not calculate the derivatives.

    double getCurv() { return _currentCurv; }
    void calcCurv();
    double getStretchedCurv() { return _stretchedCurv; }
    void calcStretchedCurv(double d); // Calculates the stretched mean curvature, and store the result in _stretchedCurv.
                                      // Does not calculate the derivatives.

};


#endif
