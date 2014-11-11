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
#include <vector>

#include "common.h"

///FORWARD DECLARATIONS
class Linker;

///MLinker class represents a cross-link between filaments

/*! The class describes interaction between 4 beads connected by a linker. Initial length of a linker is determinated by the 
 *  condition of zero initial stress, i.e., it calculated within the constructor at initiation. 
 *  A linker heads positions on a segment (between two consecutive beads on a filament)
 *  determined by two numbers (0 to 1) position1 and position2.
 */
class MLinker {
    
public:
    ///Main constructor
    ///@param position - position on cylinder 1 and 2, respectively
    ///@param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MLinker(double stretchConst, double position1, double position2,
            const vector<double>& coord11, const vector<double>& coord12,
            const vector<double>& coord21, const vector<double>& coord22);
    
    ///Getters for constants
    double getStretchingConstant(){return _kStretch;}
    double getEqLength(){return _eqLength;}
    
    ///Setter and getter for parent linker
    void setLinker(Linker* linker) {_pLinker = linker;}
    Linker* getLinker() {return _pLinker;}
    
private:
    double _eqLength;
    double _kStretch;
    
    Linker* _pLinker; ///< ptr to parent linker
    
};

#endif /* defined(__CytoMech__MLinker__) */
