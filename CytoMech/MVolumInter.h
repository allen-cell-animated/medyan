//
//  MVolumInter.h
//  CytoMech
//
//  Created by Konstantin Popov on 6/25/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MVolumInter__
#define __CytoMech__MVolumInter__

#include <iostream>
#include "Mcommon.h"

class VolumInter
{
public:
    VolumInter(Network* pn, Cylinder* pc1, Cylinder* pc2, double repusionConst);
    ~VolumInter();
    
    //Public methods called by Network:
    void CopmuteForce();
    void CopmuteForceAux();
    double ComputeEnergy();
    double ComputeEnergy(double);
    
private:
    
    // Energy calculation methods:
    double EnergyRepulsion(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4);
    double EnergyRepulsion(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double d );
    // Force calculation methods:
    void ForceRepulsion(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4 );
    void ForceRepulsionAux(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4);
   
    Cylinder* _pc1;
    Cylinder* _pc2;
    
    double _kRepuls;
    
    
    Network* _pNetwork;
    
};



#endif /* defined(__CytoMech__MVolumInter__) */
