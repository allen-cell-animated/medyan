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

void Triangle::printSelf() {
    
    cout << endl;
    
    cout << "Triangle: ptr = " << this << endl;
    cout << "Triangle ID = " << _ID << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    //cout << "Position = " << _position << endl; //< TODO: Implement this.
    
    cout << endl;
    
    /*
#ifdef CHEMISTRY
    cout << "Chemical composition of triangle:" << endl;
    _cTriangle->printCTriangle();
#endif
    */
    
    cout << endl;
    
    cout << "Vertex information..." << endl;
    
    _v[0]->printSelf();
    _v[1]->printSelf();
    _v[2]->printSelf();
    
    cout << endl;
}

