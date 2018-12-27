#ifndef MEDYAN_GTriangle_h
#define MEDYAN_GTriangle_h

#include <array>
#include <vector>

#include "MathFunctions.h"

// Forward declaration
class Triangle;

/******************************************************************************
Storing the geometric properties of the triangle patches.
******************************************************************************/

struct GTriangle {

    double area; // Current area
    double sArea; // Temporarily store the stretched area

    std::array<double, 3> _sinTheta;
    std::array<std::array<std::array<double, 3>, 3>, 3> _dSinTheta;
    std::array<double, 3> _cotTheta;
    std::array<std::array<std::array<double, 3>, 3>, 3> _dCotTheta;
    // The following variables temporarily store the angles under stretched conditions.
    std::array<double, 3> _stretchedTheta;
    std::array<double, 3> _stretchedSinTheta;
    std::array<double, 3> _stretchedCotTheta;

    mathfunc::Vec3 unitNormal; // The unit normal vector pointing outward (since the meshwork is orientable)
    mathfunc::Vec3 sUnitNormal; // Temporarily stores unit normal under stretched conditions.

    double coneVolume; // Volume of the tetrahedral formed by this triangle and the origin (0, 0, 0)
    std::array<mathfunc::Vec3, 3> dConeVolume; // The derivative of the cone volume.
    double sConeVolume;


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
    
};


#endif
