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

#include "common.h"

#include "Composite.h"
#include "CLinker.h"
#include "MLinker.h"
#include "Movable.h"
#include "Reactable.h"

///FORWARD DECLARATIONS
class Cylinder;

///Linker class is a wrapper for a chemical and mechanical linker
/*!
 * Linker class is used to create a chemical and mechanical linker when needed.
 * It contains a constructor as well as getters for mlinker and clinker.
 */

class Linker : public Composite, public Movable, public Reactable {

private:
    unique_ptr<MLinker> _mLinker; ///< ptr to mLinker
    unique_ptr<CLinker> _cLinker; ///< ptr to cLinker
    
    Cylinder* _c1; ///< first cylinder the linker is bound to
    Cylinder* _c2; ///< second cylinder the linker is bound to
    
    double _position1; ///< position on first cylinder
    double _position2; ///< position on second cylinder
    
    short _linkerType; ///integer specifying the type of linker
    int _linkerID; ///integer ID of this specific linker
    
    float _birthTime; ///Birth time of this linker
    
    Compartment* _compartment; ///< Compartment that this linker is in
    
public:
    vector<double> coordinate; ///< coordinate of midpoint, updated with updatePosition()
    
    Linker(Cylinder* c1, Cylinder* c2, short linkerType, int linkerID, double position1, double position2, bool creation);
    ~Linker();
    
    ///get cylinders
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    
    ///setter for mlinkers and clinkers
    void setCLinker(CLinker* cLinker) {_cLinker = unique_ptr<CLinker>(cLinker);}
    CLinker* getCLinker() {return _cLinker.get();}
    
    MLinker* getMLinker() {return _mLinker.get();}
    
    ///Getters and setters for position
    double getFirstPosition() {return _position1;}
    void setFirstPosition(double position1) {_position1 = position1;}
    
    double getSecondPosition() {return _position2;}
    void setSecondPosition(double position2) {_position2 = position2;}
    
    short getLinkerType() {return _linkerType;}
    int getLinkerID() {return _linkerID;}
    
    ///Update the position of this Linker
    ///@note - changes compartment of clinker if needed
    virtual void updatePosition();
    
    ///Update the reaction rates of this linker
    virtual void updateReactionRates();
    
};


#endif /* defined(__Cyto__Linker__) */
