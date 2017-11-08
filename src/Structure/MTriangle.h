#ifndef MEDYAN_MTriangle_h
#define MEDYAN_MTriangle_h

class Triangle; // Forward declaration

/*

    Storing some mechanical properties of the triangle patches

*/

class MTriangle {

private:
    Triangle* _pTriangle; // Parent triangle

    double _eqArea; // Length of unstretched area

    double _currentArea; // Current area

public:
    MTriangle(double eqArea);

    void setEqArea(double eqArea) { _eqArea = eqArea; }
    double getEqArea() { return _eqArea; }

    void setArea(double area) { _currentArea = area; }
    double getArea() { return _currentArea; }

};


#endif
