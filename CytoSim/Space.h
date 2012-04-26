//
//  SpaceShape.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_SpaceShape_h
#define CytoSim_Experimenting_SpaceShape_h

#include "Compartment.h"
#include <tuple>

class Space {
    std::vector<Compartment*> _compartments;
public:
    Space() = default;
    virtual void initializeCompartments() = 0;
    virtual void setOptions(const std::tuple<int, int, int> &options) = 0;
};

class Space1D : public Space {
private:
    int _N_Compartments;
public:
    Space1D() = default;
    void initializeCompartments();
    void setOptions(const std::tuple<int, int, int> &options) {_N_Compartments=std::get<0>(options);}
};

#endif
