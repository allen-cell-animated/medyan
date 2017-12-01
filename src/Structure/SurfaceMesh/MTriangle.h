#ifndef MEDYAN_MTriangle_h
#define MEDYAN_MTriangle_h

#include <array>

class Triangle; // Forward declaration

/*

    Storing some mechanical properties of the triangle patches

*/

class MTriangle {

private:
    Triangle* _pTriangle; // Parent triangle

    double _eqArea; // Length of unstretched area
    double _kElastic; // Elastic modulus of the triangle

    double _currentArea; // Current area
    std::array<std::array<double, 3>, 3> _dCurrentArea; // The derivative of the area. _dCurrentArea[vtxIdx][xyz]

    std::array<double, 3> _theta; // The angles corresponding to each vertex
    std::array<std::array<std::array<double, 3>, 3>, 3> _dTheta; // The derivative of the angles. _dTheta[angleIdx][vtxIdx][xyz]
                                                          // For example, _dTheta[0][2][1] means d(_theta[0]) / dy for vertex 2.
    std::array<double, 3> _sinTheta;
    std::array<std::array<std::array<double, 3>, 3>, 3> _dSinTheta;
    std::array<double, 3> _cotTheta;
    std::array<std::array<std::array<double, 3>, 3>, 3> _dCotTheta;

public:
    MTriangle(double eqArea);

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setElasticModulus(double kElastic) { _kElastic = kElastic; }
    double getElasticModulus() { return _kElastic; }

    void setArea(double area) { _currentArea = area; }
    double getArea() { return _currentArea; }

    std::array<double, 3>& getTheta() { return _theta; }
    void calcTheta();

};


#endif
