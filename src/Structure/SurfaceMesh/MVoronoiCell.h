#ifndef MEDYAN_MVoronoiCell_h
#define MEDYAN_MVoronoiCell_h

#include <vector>

class Vertex; // Forward declaration

class MVoronoiCell {

private:
    Vertex* _vertex; // Parent triangle

    double _eqArea; // Length of unstretched area, determined at mesh generation
    double _kElastic; // Elastic modulus of the triangle

    double _currentArea; // Current area


public:
    MVoronoiCell(double eqArea) {
        setEqArea(eqArea);
    }

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setElasticModulus(double kElastic) { _kElastic = kElastic; }
    double getElasticModulus() { return _kElastic; }

    void setArea(double area) { _currentArea = area; }
    double getArea() { return _currentArea; }

};


#endif
