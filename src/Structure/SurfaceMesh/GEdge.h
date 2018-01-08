#ifndef MEDYAN_GEdge_h
#define MEDYAN_GEdge_h

#include <array>

class Edge; // Forward declaration

class GEdge {

private:
    Edge* _pEdge; // Parent edge

    double _currentLength; // Current length
    std::array<std::array<double, 3>, 2> _dCurrentLength; // The derivative of the length. _dCurrentLength[vtxIdx][xyz]
    double _stretchedLength; // Temporarily store the stretched length

    std::array<double, 3> _pseudoUnitNormal; // The pseudo unit normal vector at the edge pointing outward.
    std::array<double, 3> _stretchedPseudoUnitNormal; // The pseudo normal under stretch

public:
    GEdge() {}

    /// Set parent 
    void setEdge(Edge* e) { _pEdge = e; }
    Edge* getEdge() { return _pEdge; }

    double getLength() { return _currentLength; }
    std::array<std::array<double, 3>, 2>& getDLength() { return _dCurrentLength; }
    void calcLength();
    double getStretchedLength() { return _stretchedLength; }
    void calcStretchedLength(double d); // Calculates the stretched length, and store the result in _stretchedLength.
                                        // Does not calculate the derivatives.
    
    std::array<double, 3>& getPseudoUnitNormal() { return _pseudoUnitNormal; }
    void calcPseudoUnitNormal(); // Calculates the pseudo unit normal w/o derivatives.
    std::array<double, 3>& getStretchedPseudoUnitNormal() { return _stretchedPseudoUnitNormal; }
    void calcStretchedPseudoUnitNormal(double d); // Calculates the pseudo unit normal w/o derivatives.
    
};


#endif
