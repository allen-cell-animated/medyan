
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_LinkerStretching_h
#define MEDYAN_LinkerStretching_h

#include "common.h"

#include "LinkerInteractions.h"
#include "Linker.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents a Linker stretching interaction
template <class LStretchingInteractionType>
class LinkerStretching : public LinkerInteractions {
    
private:
    LStretchingInteractionType _FFType;

    int *beadSet;
    
    ///Array describing the constants in calculation
    double *kstr;
    double *eql;
    double *pos1;
    double *pos2;
    
public:
    
    ///Array describing indexed set of interactions
    ///For linkers, this is a 4-bead potential
    const static int n = 4;
    
    ///< Constructor initializes data
    LinkerStretching () {
        beadSet = new int[n * Linker::getLinkers().size()];
        kstr = new double[Linker::getLinkers().size()];
        eql = new double[Linker::getLinkers().size()];
        pos1 = new double[Linker::getLinkers().size()];
        pos2 = new double[Linker::getLinkers().size()];
    }
    
    ~LinkerStretching () { delete beadSet; delete kstr;
                           delete eql; delete pos1; delete pos2;}
    
    virtual void vectorizeInteractions();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual const string getName() {return "Linker Stretching";}
};

#endif
