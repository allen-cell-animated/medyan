//
//  System.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_System_h
#define CytoSim_Experimenting_System_h

#include <vector>
#include "Reaction.h"
#include "Species.h"
#include "Compartment.h"
#include "Space.h"
#include "ParseSystem.h"

class System
{
private:
    std::vector<Compartment*> _compartments;
    CompartmentProto _comp_proto;
    Space* _space;
    ParseSystem _parse_system;
    //Note: bulk species and their reactions should go here
public:
    System() : _comp_proto(this) {
    }
    void setSpace(Space* space){_space=space;}
    void setSpaceOptions(const std::tuple<int, int, int> &options){_space->setOptions(options);}
    void initializeCompartments() {_space->initializeCompartments();}
    bool parseSpeciesTypes(){
        _parse_system.parseSpecies();
        return true;
    }
    bool parseSystem(){
        _comp_proto.setSpecies();
        _comp_proto.printSpecies();
        return true;
    }
    //Accessors

    SpeciesType* getSpeciesType(const std::string &s){return _parse_system.getSpeciesType(s);}
    std::vector<SpeciesType*> allSpeciesType(){return _parse_system.allSpeciesType();}
    int getNumCompartments() const {return _compartments.size();}

};

#endif
