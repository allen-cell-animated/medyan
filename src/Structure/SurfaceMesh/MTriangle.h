#ifndef MEDYAN_MTriangle_h
#define MEDYAN_MTriangle_h

#include <array>

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

class MTriangle {

private:

    double _kExVol; ///< Local excluded volume constant, which describes
                    ///< excluded volume interactions between triangles and cylinder tips

public:
    MTriangle(short membraneType);

    void getExVolConst(double kExVol) { _kExVol = kExVol; }
    double getExVolConst() { return _kExVol; }

};


#endif
