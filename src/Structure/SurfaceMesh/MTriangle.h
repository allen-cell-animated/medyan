#ifndef MEDYAN_MTriangle_h
#define MEDYAN_MTriangle_h

#include <array>

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

class MTriangle {

private:

    [[deprecated]] double _eqArea; // Length of unstretched area
    [[deprecated]] double _kElastic; // Elastic modulus of the triangle

    double _kExVol; ///< Local excluded volume constant, which describes
                    ///< excluded volume interactions between triangles and cylinder tips

public:
    MTriangle(short membraneType);

    void getExVolConst(double kExVol) { _kExVol = kExVol; }
    double getExVolConst() { return _kExVol; }

};


#endif
