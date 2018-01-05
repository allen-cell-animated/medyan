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

public:
    MVoronoiCell(short membraneType);

    /// Set parent 
    void setVertex(Vertex* v) { _pVertex = v; }
    Vertex* getVertex() { return _pVertex; }

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setElasticModulus(double kElastic) { _kElastic = kElastic; }
    double getElasticModulus() { return _kElastic; }

    void setBendingModulus(double kBending) { _kBending = kBending; }
    double getBendingModulus() { return _kBending; }

    void setEqCurv(double eqCurv) { _eqCurv = eqCurv; }
    double getEqCurv() { return _eqCurv; }

};


#endif
