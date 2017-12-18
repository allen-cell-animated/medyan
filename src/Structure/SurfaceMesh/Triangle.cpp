#include "Triangle.h"

Database<Triangle*> Triangle::_triangles;

Triangle::Triangle(Composite* parent, Vertex* v1, Vertex* v2, Vertex* v3):
    Trackable(true, false, true, false),
    _v{v1, v2, v3}, _edges{nullptr, nullptr, nullptr}, _id(_triangles.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mTriangle = unique_ptr<MTriangle>(new MTriangle(getType()));
    _mTriangle->setTriangle(this);
#endif

}

void Triangle::updateCoordinate() {
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (_v[0]->coordinate[coordIdx] + _v[1]->coordinate[coordIdx] + _v[2]->coordinate[coordIdx]) / 3;
    }
}

void Triangle::updatePosition() {
    updateCoordinate();
    
    // If compartments are implemented, they are also updated here.
}

void Triangle::printSelf() {
    
    cout << endl;
    
    cout << "Triangle: ptr = " << this << endl;
    cout << "Triangle ID = " << _id << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    //cout << "Position = " << _position << endl; //< TODO: Implement this.
    
    cout << endl;
    
    /*
#ifdef CHEMISTRY
    cout << "Chemical composition of triangle:" << endl;
    _cTriangle->printCTriangle();
#endif
    
    cout << endl;
    */
    
    cout << "Vertex information..." << endl;
    
    _v[0]->printSelf();
    _v[1]->printSelf();
    _v[2]->printSelf();
    
    cout << endl;
}

void Triangle::updateGeometry(bool calcDerivative, double d) {

#ifdef MECHANICS
    _mTriangle->updateGeometry(calcDerivative, d);
#endif

}
