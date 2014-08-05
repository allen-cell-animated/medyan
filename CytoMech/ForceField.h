//
//  ForceField.h
//  Cyto
//
//  Created by James Komianos on 8/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ForceField__
#define __Cyto__ForceField__

#include <iostream>

///ForceField is an abstract class to represent various force field calculations
/*!
 *  ForceField is used for force calculations between filaments, beads, linkers, etc.
 *  Specific implementations of the ForceField class will have different potentials.
 */
class ForceField {

private:
    std::string _name;
    
public:
    ///Constructor
    ForceField(std::string name) {_name = name;}
    ~ForceField();
    
    std::string getName() {return _name;}
    
    virtual double ComputeEnergy() = 0;
    
    virtual void ComputeForces() = 0;
    
};






#endif /* defined(__Cyto__ForceField__) */
