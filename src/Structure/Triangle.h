#ifndef MEDYAN_Triangle_h
#define MEDYAN_Triangle_h

#include <array>

#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Bead.h"

class Triangle:
    public Component,
    public Trackable,
    public Movable,
    public Reactable,
    public DynamicNeighbor {

private:
    // Pointers to 3 beads. The rotation of beads must be pointing outside.
    std::array<Bead*, 3> _b;

public:
    Triangle(Composite *parent, Bead *b1, Bead *b2, Bead *b3);

};


#endif
