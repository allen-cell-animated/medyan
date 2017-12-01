#ifndef MEDYAN_MEdge_h
#define MEDYAN_MEdge_h

#include <array>

class Edge; // Forward declaration

class MEdge {

private:
    Edge* _pEdge; // Parent edge

    double _currentLength; // Current length
    std::array<std::array<double, 3>, 2> _dCurrentLength; // The derivative of the length. _dCurrentLength[vtxIdx][xyz]

public:
    MEdge() {}

    double getLength() { return _currentLength; }
    void calcLength();
};


#endif
