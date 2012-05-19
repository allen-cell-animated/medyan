//
//  ChemSim.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemSim_h
#define CytoSim_ChemSim_h

#include <memory>
#include "Reaction.h"

namespace chem {

class ChemSimImpl {
public:
    virtual ~ChemSimImpl() {};
    virtual void initialize() = 0;
    virtual void addReaction(Reaction *r) = 0;
    virtual void removeReaction(Reaction *r) = 0;
    virtual void run(int steps) = 0;
    virtual void printReactions() const = 0;
};

class ChemSim {
public:
    // Constructor - no ownership of ChemSimImpl ptr
    ChemSim(ChemSimImpl *csi);
    ChemSim(const ChemSim &rhs) = delete;
    ChemSim& operator=(ChemSim &rhs) = delete;
    void initialize();
    void addReaction(Reaction *r);
    void removeReaction(Reaction *r);
    void run(int steps);
    void printReactions() const;
private:
    ChemSimImpl* _pimpl;
};

} // end of namespace
#endif
