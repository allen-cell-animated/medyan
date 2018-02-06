#ifndef MEDYAN_GMembrane_h
#define MEDYAN_GMembrane_h

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
Stores geometric properties of the membrane.
******************************************************************************/

class GMembrane {

private:

    Membrane* _pMembrane; // Parent membrane

    double _volume; // The volume of the enclosed membrane

public:
    // Constructors
    GMembrane() {}

    // Getters and setters
    void setMembrane(Membrane* m) { _pMembrane = m; }
    Membrane* getMembrane()const { return _pMembrane; }

    double getVolume()const { return _volume; }
    void calcVolume();
    
};





#endif
