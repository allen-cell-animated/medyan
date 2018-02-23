#ifndef MEDYAN_MMembrane_h
#define MEDYAN_MMembrane_h

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

class MMembrane {

private:

    Membrane* _pMembrane; // Parent membrane

    double _eqVolume; // The equilibrium volume

public:
    // Constructors
    MMembrane() {};

    // Getters and setters
    void setMembrane(Membrane* m) { _pMembrane = m; }
    Membrane* getMembrane()const { return _pMembrane; }

    double getEqVolume()const { return _eqVolume; }
    void setEqVolume(double eqVolume) { _eqVolume = eqVolume; }
    
};





#endif
