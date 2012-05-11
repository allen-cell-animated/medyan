//
//  ChemNRM.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemNRM_h
#define CytoSim_ChemNRM_h

#include <memory>
#include "Reaction.h"

class ChemNRMImpl;

class ChemNRM {
public:
    ChemNRM();
    ChemNRM(const ChemNRM &rhs) = delete;
    ChemNRM& operator=(ChemNRM &rhs) = delete;
    void initialize();
    void addReaction(Reaction *r);
    void removeReaction(Reaction *r);
    void run(int steps);
    size_t getSize() const;
    float getTime() const;
    void printReactions() const;
private:
    std::unique_ptr<ChemNRMImpl> _pimpl;
};

#endif
