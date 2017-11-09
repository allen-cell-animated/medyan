#ifndef MEDYAN_Membrane_h
#define MEDYAN_Membrane_h

#include <vector>

#include "Trackable.h"
#include "Composite.h"

// FORWARD DECLARATIONS
class SubSystem;
class Triangle;
class Bead;

class Membrane: public Composite, public Trackable {

    friend class Controller;

private:

    SubSystem* _subSystem; // SubSystem pointer
};








#endif
