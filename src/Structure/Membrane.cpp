#include "Membrane.h"

#include "common.h"

#include "SubSystem.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(SubSystem* s, short membraneType,
    vector<tuple<array<double, 3>, vector<size_t>>>& membraneData):
    Trackable(), _subSystem(s), _memType(membraneType), _Id(_membranes.getID()) {
    
    // Build the meshwork using vertex and neighbor information

    // TODO: Implement this

}

Membrane::~Membrane() {
    for(auto& v: _vertexVector) _subSystem->removeTrackable<Vertex>(v);
    for(auto& e: _edgeVector) _subSystem->removeTrackable<Edge>(e);
    for(auto& t: _triangleVector) _subSystem->removeTrackable<Triangle>(t);
}

void Membrane::printSelf() {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    cout << "Membrane Id = " << _Id << endl;
    cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    cout << "Triangle information..." << endl;
    
    for(auto t : _triangleVector)
        t->printSelf();
    
    cout << endl;
    
}
