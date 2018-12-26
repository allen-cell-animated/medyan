#ifndef MEDYAN_GMembrane_h
#define MEDYAN_GMembrane_h

#include <array>
#include <vector>

#include "MathFunctions.h"

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

struct GMembrane {

    double volume; // The volume of the enclosed membrane
    std::vector<mathfunc::Vec3> dVolume; // The derivative of volume on vertices.
    double sVolume; // The stretched volume

    // Constructors
    GMembrane() {}

    // Getters and setters
    void setMembrane(Membrane* m) { _pMembrane = m; }
    Membrane* getMembrane()const { return _pMembrane; }

    double getVolume()const { return _volume; }
    std::vector<std::array<double, 3>>& getDVolume() { return _dVolume; }
    void calcVolume();
    double getStretchedVolume()const { return _stretchedVolume; }
    void calcStretchedVolume(double d); // Does not calculate derivative
    
};





#endif
