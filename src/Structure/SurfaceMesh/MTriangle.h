#ifndef MEDYAN_MTriangle_h
#define MEDYAN_MTriangle_h

#include <array>

class Triangle; // Forward declaration

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

class MTriangle {

private:
    Triangle* _pTriangle; // Parent triangle

    double _eqArea; // Length of unstretched area
    double _kElastic; // Elastic modulus of the triangle

    double _kExVol; ///< Local excluded volume constant, which describes
                    ///< excluded volume interactions between triangles and cylinder tips

    double _currentArea; // Current area
    std::array<std::array<double, 3>, 3> _dCurrentArea; // The derivative of the area. _dCurrentArea[vtxIdx][xyz]
    double _stretchedArea; // Temporarily store the stretched area

    std::array<double, 3> _theta; // The angles corresponding to each vertex
    std::array<std::array<std::array<double, 3>, 3>, 3> _dTheta; // The derivative of the angles. _dTheta[angleIdx][vtxIdx][xyz]
                                                          // For example, _dTheta[0][2][1] means d(_theta[0]) / dy for vertex 2.
    std::array<double, 3> _sinTheta;
    std::array<std::array<std::array<double, 3>, 3>, 3> _dSinTheta;
    std::array<double, 3> _cotTheta;
    std::array<std::array<std::array<double, 3>, 3>, 3> _dCotTheta;
    // The following variables temporarily store the angles under stretched conditions.
    std::array<double, 3> _stretchedTheta;
    std::array<double, 3> _stretchedSinTheta;
    std::array<double, 3> _stretchedCotTheta;

public:
    MTriangle(short membraneType);

    /// Set parent 
    void setTriangle(Triangle* t) { _pTriangle = t; }
    Triangle* getTriangle() { return _pTriangle; }

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setElasticModulus(double kElastic) { _kElastic = kElastic; }
    double getElasticModulus() { return _kElastic; }

    void getExVolConst(double kExVol) { _kExVol = kExVol; }
    double getExVolConst() { return _kExVol; }

    // Not allowing setting the area: void setArea(double area) { _currentArea = area; }
    double getArea() { return _currentArea; }
    std::array<std::array<double, 3>, 3>& getDArea() { return _dCurrentArea; }
    void calcArea();
    double getStretchedArea() { return _stretchedArea; }
    void calcStretchedArea(double d); // Calculates the stretched area, and store the result in _stretchedArea
                                      // Does not calculate derivatives.

    std::array<double, 3>& getTheta() { return _theta; }
    std::array<double, 3>& getSinTheta() { return _sinTheta; }
    std::array<double, 3>& getCotTheta() { return _cotTheta; }
    std::array<std::array<std::array<double, 3>, 3>, 3>& getDTheta() { return _dTheta; }
    std::array<std::array<std::array<double, 3>, 3>, 3>& getDSinTheta() { return _dSinTheta; }
    std::array<std::array<std::array<double, 3>, 3>, 3>& getDCotTheta() { return _dCotTheta; }
    void calcTheta(); // would calculate all variables related to theta
    std::array<double, 3>& getStretchedTheta() { return _stretchedTheta; }
    std::array<double, 3>& getStretchedSinTheta() { return _stretchedSinTheta; }
    std::array<double, 3>& getStretchedCotTheta() { return _stretchedCotTheta; }
    void calcStretchedTheta(double d); // Calculates angles under stretched conditions (w/o derivatives)
                                       // The results are stored in _stretchedXxx variables.

};


#endif
