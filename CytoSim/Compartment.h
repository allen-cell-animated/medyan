//
//  Compartment.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Compartment_h
#define CytoSim_Experimenting_Compartment_h

#include <vector>
#include "Reaction.h"
#include "Species.h"

class System;

class Compartment {
private:
    std::vector<Species*> _species;
    std::vector<ReactionBase*> _reactions;
    compartment_num_t _id;
public:
    Compartment()=default;
};

class CompartmentProto {
private:
    //Compartment _comp;
    std::vector<Species*> _species_proto;
    std::vector<ReactionBase*> _reactions_proto;
    System *_system;
public:
    CompartmentProto(System *system) : _system(system) {};
    bool setSpecies();
    void printSpecies(){
        for(auto s: _species_proto)
            s->printSelf();
    }
};

#endif
