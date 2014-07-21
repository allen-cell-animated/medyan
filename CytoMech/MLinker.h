//
//  MLinker.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MLinker__
#define __CytoMech__MLinker__

#include <iostream>
#include "Mcommon.h"



class Linker{
    
    /// The class which manages interactions betwin cross-linkered beads.
    
public:
    Linker(Network* pn, Cylinder* pc1, Cylinder* pc2, double stretchConst);
    //Public methods called by Network:
    void CopmuteForce();
    void CopmuteForceAux();
    double CopmuteEnergy();
    double CopmuteEnergy(double);
    
//    //More for test purpose:
//    Bead* getFirstBead(){return _pb1;}
//    Bead* getSecondBead(){return _pb2;}
    
    
    
private:
    
    // Energy calculation methods:
    double EnergyHarmonicStretching(Bead* pb1, Bead* pb2 );
    double EnergyHarmonicStretching(Bead* pb1, Bead* pb2, double d );
    // Force calculation methods:
    void ForceHarmonicStretching(Bead* pb1, Bead* pb2 );
    void ForceHarmonicStretchingAux(Bead* pb1, Bead* pb2);
    
    Cylinder* _pc1;
    Cylinder* _pc2;
    double _eqLength;
    double _kStretch;
    
    
    Network* _pNetwork;
    
};

#endif /* defined(__CytoMech__MLinker__) */
