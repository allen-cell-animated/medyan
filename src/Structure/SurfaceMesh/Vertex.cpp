#include "Vertex.h"

Database<Vertex*> Vertex::_vertices;

Vertex::Vertex(vector<double> v, Composite* parent, int position):
    Bead(v, parent, position) {}

Vertex::Vertex(Composite* parent, int position):
    Bead(parent, position) {}
