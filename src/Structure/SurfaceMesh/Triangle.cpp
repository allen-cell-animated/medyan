#include "Triangle.h"

Database<Triangle*> Triangle::_triangles;

Triangle::Triangle(Composite* parent, Vertex* v1, Vertex* v2, Vertex* v3):
    _v{v1, v2, v3}, _edges{nullptr, nullptr, nullptr} {
    
    parent -> addChild(unique_ptr<Component>(this));

}

void Triangle::updateCoordinate() {
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (_v[0]->coordinate[coordIdx] + _v[1]->coordinate[coordIdx] + _v[2]->coordinate[coordIdx]) / 3;
    }
}
