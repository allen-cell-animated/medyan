#ifndef MEDYAN_GTriangle_h
#define MEDYAN_GTriangle_h

#include <array>

class Triangle; // Forward declaration

/******************************************************************************

Storing the geometric properties of the triangle patches

******************************************************************************/

class GTriangle {

private:
    Triangle* _pTriangle; // Parent triangle

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

    std::array<double, 3> _unitNormal; // The unit normal vector pointing outward (since the meshwork is orientable)
    std::array<double, 3> _stretchedUnitNormal; // Temporarily stores unit normal under stretched conditions.

public:
    GTriangle() {}

    /// Set parent 
    void setTriangle(Triangle* t) { _pTriangle = t; }
    Triangle* getTriangle() { return _pTriangle; }

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
    
    std::array<double, 3>& getUnitNormal() { return _unitNormal; }
    void calcUnitNormal(); // Calculates the unit normal of the triangle (w/o derivatives)
    std::array<double, 3>& getStretchedUnitNormal() { return _stretchedUnitNormal; }
    void calcStretchedUnitNormal(double d); // Calculates the unit normal under stretched conditions (w/o derivatives)
    
    void updateGeometry(bool calcDerivative=false, double d=0.0) {
        // Currently, derivative cannot be calculated for d != 0
        if(calcDerivative) {
            calcTheta();
            calcArea();
            calcUnitNormal();
        }
        else {
            calcStretchedTheta(d);
            calcStretchedArea(d);
            calcStretchedUnitNormal(d);
        }
    }

};


#endif
