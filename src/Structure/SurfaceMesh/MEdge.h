#ifndef MEDYAN_MEdge_h
#define MEDYAN_MEdge_h

#include <array>

class Edge; // Forward declaration

class MEdge {

private:
    Edge* _pEdge; // Parent edge

    double _currentLength; // Current length
    std::array<std::array<double, 3>, 2> _dCurrentLength; // The derivative of the length. _dCurrentLength[vtxIdx][xyz]
    double _stretchedLength; // Temporarily store the stretched length

public:
    MEdge() {}

    /// Set parent 
    void setEdge(Edge* e) { _pEdge = e; }
    Edge* getEdge() { return _pEdge; }

    double getLength() { return _currentLength; }
    std::array<std::array<double, 3>, 2>& getDLength() { return _dCurrentLength; }
    void calcLength();
    double getStretchedLength() { return _stretchedLength; }
    void calcStretchedLength(double d); // Calculates the stretched length, and store the result in _stretchedLength.
                                        // Does not calculate the derivatives.
    
    void updateGeometry(bool calcDerivative=false, double d=0.0) {
        // Currently, derivative cannot be calculated for d != 0
        if(calcDerivative) calcLength();
        else calcStretchedLength(d);
    }
};


#endif
