#ifndef MEDYAN_SimpleGeometricObject_h
#define MEDYAN_SimpleGeometricObject_h

#include <array>
#include <stdexcept>
#include <vector>

#include "Geometric.h"
#include "MathFunctions.h"
using namespace mathfunc;

// Forward declarations

/******************************************************************************

Simple geometric objects and some math functions. As the name suggests, these
objects do not participate in force calculations.

******************************************************************************/

//The simple vertex is a point of arbitrary dimension.
template<size_t Dim>
class SimpleVertex: public Geometric {
private:
    std::array<double, Dim> _coordinate;

public:

    // Main constructor
    SimpleVertex(const std::array<double, Dim>& newCoordinate): _coordinate(newCoordinate) {}

    // Properties
    const std::array<double, Dim>& getCoordinate()const { return _coordinate; }
    void setCoordinate(const std::array<double, Dim>& newCoordinate) {
        _coordinate = newCoordinate;
    }
    double operator[](size_t n)const { return _coordinate[n]; }
    double& operator[](size_t n) { return _coordinate[n]; }

    // Geometry
    virtual void updateGeometry()override {}
};

inline SimpleVertex<3> changeDimension_2_3(const SimpleVertex<2>& v, size_t aspect, double value) {
    // aspect 0: yOz, 1: zOx, 2: xOy
    switch(aspect) {
        case 0: return SimpleVertex<3>{{value, v[0], v[1]}};
        case 1: return SimpleVertex<3>{{v[1], value, v[0]}};
        case 2: return SimpleVertex<3>{{v[0], v[1], value}};
        default: throw out_of_range("Aspect out of range.");
    }
}
inline SimpleVertex<2> changeDimension_3_2(const SimpleVertex<3>& v, size_t aspect) {
    // aspect 0: yOz, 1: zOx, 2: xOy
    switch(aspect) {
        case 0: return SimpleVertex<2>{{v[1], v[2]}};
        case 1: return SimpleVertex<2>{{v[2], v[0]}};
        case 2: return SimpleVertex<2>{{v[0], v[1]}};
        default: throw out_of_range("Aspect out of range.");
    }
}

/******************************************************************************
The simple polygon is made of closed connected vertices.

Assumptions (not checked):
    - all the vertices are on the same plane
    - the polygon does not intersect itself, implying orientability and no hole

For 2D polygon, vertices are connected in the counter-clockwise direction
******************************************************************************/
template<size_t Dim>
class SimplePolygon: public Geometric {
private:
    double _area;

    std::vector<SimpleVertex<Dim>*> _vertices;

public:
    // Main constructor
    SimplePolygon();
    SimplePolygon(const std::vector<SimpleVertex<Dim>*>& newVertices): _vertices(newVertices) {}

    // Properties
    double getArea()const { return _area; }
    virtual void calcArea();

    const SimpleVertex<Dim>& vertex(size_t n)const { return *(_vertices[n]); }

    void addVertex(const SimpleVertex<Dim>* v, size_t loc=0) {
        _vertices.insert(_vertices.begin() + loc, v);
    }

    // Geometry
    virtual void updateGeometry()override {
        calcArea();
    }
}

template<>
inline double SimplePolygon<3>::calcArea() {
    size_t n = vertices.size();
    if(n <= 2) { _area = 0; return; }

    auto normal = normalizedVector(vectorProduct(
        vertex(0).getCoordinate(), vertex(1).getCoordinate(),
        vertex(0).getCoordinate(), vertex(2).getCoordinate()
    ));

    std::array<double, 3> tempCross {};
    for(size_t idx = 0; idx < n; ++idx) {
        vectorIncrease(tempCross, crossProduct(vertex(idx).getCoordinate(), vertex((idx + 1) % n).getCoordinate()));
    }

    _area = 0.5 * dotProduct(normal, tempCross);
    
}
template<>
inline double SimplePolygon<2>::calcArea() {
    size_t n = vertices.size();
    if(n <= 2) { _area = 0; return; }

    _area = 0;
    for(size_t idx = 0; idx < n; ++idx) {
        _area += vertex(idx)[0] * vertex((idx + 1) % n)[1] - vertex((idx + 1) % n)[0] * vertex(idx)[1];
    }
    _area /= 2;

}



#endif
