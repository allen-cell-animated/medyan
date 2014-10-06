//
//  Linker.h
//  Cyto
//
//  Created by James Komianos on 10/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Linker__
#define __Cyto__Linker__

#include <iostream>
#include "Composite.h"
#include "CLinker.h"
#include "MLinker.h"
#include "Cylinder.h"

class Compartment;
class SpeciesBound;

///Linker class is a wrapper for a chemical and mechanical linker
/*!
 * Linker class is used to create a chemical and mechanical linker when needed.
 * It contains a constructor as well as getters for mlinker and clinker.
 */

class Linker : public Composite {

private:
    std::unique_ptr<MLinker> _mLinker; ///< ptr to mLinker
    std::unique_ptr<CLinker> _cLinker; ///< ptr to cLinker
    
    Cylinder* _pc1; ///< first cylinder the linker is bound to
    Cylinder* _pc2; ///< second cylinder the linker is bound to
    
    double _position1; ///< position on first cylinder
    double _position2; ///< position on second cylinder
    
public:
    Linker(Cylinder* pc1, Cylinder* pc2, Compartment* c, double position1, double position2);
    ~Linker() {}
    
    ///get cylinders
    Cylinder* getFirstCylinder() {return _pc1;}
    Cylinder* getSecondCylinder() {return _pc2;}
    
    ///setter for mlinkers and clinkers
    void setCLinker(CLinker* cLinker) {_cLinker = std::unique_ptr<CLinker>(cLinker);}
    CLinker* getCLinker() {return _cLinker.get();}
    
    void setMLinker(MLinker* mLinker) {_mLinker = std::unique_ptr<MLinker>(mLinker);}
    MLinker* getMLinker() {return _mLinker.get();}
    
    ///Getters and setters for position
    double getFirstPosition() {return _position1;}
    void setFirstPosition(double position1) {_position1 = position1;}
    double getSecondPosition() {return _position2;}
    void setSecondPosition(double position2) {_position2 = position2;}
    
    ///Update the position of this Linker
    ///@note - changes compartment of clinker if needed
    void updatePosition();
    
};


#endif /* defined(__Cyto__Linker__) */
