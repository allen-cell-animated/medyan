#ifndef MEDYAN_MMembrane_h
#define MEDYAN_MMembrane_h

// FORWARD DECLARATIONS
class Membrane;

/******************************************************************************
Stores mechanical properties of the membrane.
******************************************************************************/

class MMembrane {

private:

    Membrane* _pMembrane; // Parent membrane

    double _kComp; // Compressibility constant

    double _volume; // The volume of the enclosed membrane

public:
    // Constructors
    Membrane();

    // Getters and setters
    void setMembrane(Membrane* m) { _pMembrane = m; }
    Membrane* getMembrane()const { return _pMembrane; }
    
    void setCompConst(double kComp) { _kComp = kComp; }
    double getCompConst()const { return _kComp; }

    double getVolume()const { return _volume; }
    void calcVolume();
    
};





#endif
