#ifndef MEDYAN_MTriangle_h
#define MEDYAN_MTriangle_h

#include <array>

class Triangle; // Forward declaration

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

class MTriangle {

private:
    Triangle* _pTriangle; // Parent triangle

    double _eqArea; // Length of unstretched area
    double _kElastic; // Elastic modulus of the triangle

    double _kExVol; ///< Local excluded volume constant, which describes
                    ///< excluded volume interactions between triangles and cylinder tips

public:
    MTriangle(short membraneType);

    /// Set parent 
    void setTriangle(Triangle* t) { _pTriangle = t; }
    Triangle* getTriangle() { return _pTriangle; }

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setElasticModulus(double kElastic) { _kElastic = kElastic; }
    double getElasticModulus() { return _kElastic; }

    void getExVolConst(double kExVol) { _kExVol = kExVol; }
    double getExVolConst() { return _kExVol; }

};


#endif
