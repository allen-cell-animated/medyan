#ifndef MEDYAN_GMembrane_h
#define MEDYAN_GMembrane_h

#include <array>
#include <vector>

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

class GMembrane {

private:

    Membrane* _pMembrane; // Parent membrane

    double _volume; // The volume of the enclosed membrane
    std::vector<std::array<double, 3>> _dVolume; // The derivative of volume on vertices.
    double _stretchedVolume; // The stretched volume

public:
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
