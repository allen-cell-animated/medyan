//
//  MFilament.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MFilament__
#define __CytoMech__MFilament__

#include <iostream>
#include <iostream>
#include <deque>
#include "Mcommon.h"
#include "MComposite.h"
#include <math.h>

class CylinderDBKey;
class BeadDBKey;

class Filament : public MComposite
{
    /*! 
     *This class contein information about beads connectivity. It iterates over the beads whithin itself and compute a bonede contribution to forces and energy of beads.
     */
    
private:
    std::deque<Cylinder*> _pCylinderVector; //< Vector of cylinders;
    Cylinder* _pLastCylinder;
    System* _pSubSystem;
    Network* _pNetwork;

    
public:
    
	Filament(System* ps, Network* pn, std::vector<double> position, std::vector<double> direction);
    Filament(System* ps, Network* pn, std::vector<std::vector<double>> position, int numBeads);
    
    void PolymerizeFront();  // Polymerization of a new cylinder (with the first bead = new bead b and last bead = empty: before: --x---x---o, after --x---x---[x---o]). Next position is based on previous beads directions in the filament. This function creates a new bead. So, this function mostly called during further polymerization, not initiation.

    void PolymerizeBack();  // Same as Polymerization frond, but adds a new first cylinder with first bead = new bead and a second bead is equal to the firs bead in the cylynder, which used to be first.
    void PolymerizeFront(std::vector<double> coordonates );
    void PolymerizeBack(std::vector<double> coordonates );
    
    void DeleteBead(Bead*);
    
    Cylinder* getLastCylinder() {return _pLastCylinder;}
    
    std::deque<Cylinder*> getCylinderVector() {return _pCylinderVector;}
    
    
    
    //just for example, will rewrite this function, so it not returns anything
    std:: vector<double> NextBeadProjection(Bead* pb, double d, std::vector<double> director);
    std::vector<std::vector<double> > StraightFilamentProjection(std::vector<std::vector<double>> v, int numBeads);
    std:: vector<std::vector<double> > ArcFilamentProjection(std::vector<std::vector<double>> v, int numBeads);
};

#endif /* defined(__CytoMech__MFilament__) */
