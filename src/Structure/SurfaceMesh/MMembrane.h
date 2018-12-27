#ifndef MEDYAN_MMembrane_h
#define MEDYAN_MMembrane_h

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

class MMembrane {

private:

    double _eqVolume; // The equilibrium volume

    double _eqArea; // The equilibrium area
    double _kElastic; // Area elasticity

public:

    double getEqVolume()const { return _eqVolume; }
    void setEqVolume(double eqVolume) { _eqVolume = eqVolume; }

    double getEqArea()const { return _eqArea; }
    double getKElastic()const { return _kElastic; }
    
};





#endif
