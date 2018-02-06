#include "GMembrane.h"

#include "Membrane.h"
#include "Triangle.h"
#include "MTriangle.h"
#include "Vertex.h"

void GMembrane::calcVolume() {
}
// TODO
class MMembrane {

private:

    Membrane* _pMembrane; // Parent membrane

    double _kComp; // Compressibility constant
    double _eqVolume; // Equilibrium volume

    double _volume; // The volume of the enclosed membrane

public:
    // Constructors
    Membrane();

    // Getters and setters
    void setMembrane(Membrane* m) { _pMembrane = m; }
    Membrane* getMembrane()const { return _pMembrane; }

    void setCompConst(double kComp) { _kComp = kComp; }
    double getCompConst()const { return _kComp; }

    void setEqVolume(double eqVolume) { _eqVolume = eqVolume; }
    double getEqVolume()const { return _eqVolume; }

    double getVolume()const { return _volume; }
    void calcVolume();
    
};





#endif
