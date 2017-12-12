#include "Membrane.h"

Database<Membrane*> Membrane::_membranes;

void Membrane::printSelf() {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    //cout << "Membrane ID = " << _ID << endl;
    //cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    cout << "Triangle information..." << endl;
    
    for(auto t : _triangleVector)
        t->printSelf();
    
    cout << endl;
    
}
