#ifndef MEDYAN_MVoronoiCell_h
#define MEDYAN_MVoronoiCell_h

class MVoronoiCell {

private:

    [[deprecated]] double _eqArea; // Length of unstretched area, determined at mesh generation
    [[deprecated]] double _kElastic; // Local elastic modulus
    double _kBending; // Local bending modulus
    double _eqCurv; // Local spontaneous curvature

public:
    MVoronoiCell(short membraneType);

    void setBendingModulus(double kBending) { _kBending = kBending; }
    double getBendingModulus() { return _kBending; }

    void setEqCurv(double eqCurv) { _eqCurv = eqCurv; }
    double getEqCurv() { return _eqCurv; }

};


#endif
