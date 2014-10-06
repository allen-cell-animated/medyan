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

class Linker;

///MLinker class represents a cross-link between filaments
/*!
 *  A MLinker represents a cross link between cylinders. It contains mechanical information, including
 *  constants and a back-pointer to the parent linker object.
 */
class MLinker {
    
public:
    ///Main constructor
    ///@param position - position on cylinder 1 and 2, respectively
    ///@param coord - coordinates of cylinder1's bead 1, bead 2, etc
    MLinker(double stretchConst, double position1, double position2,
            const std::vector<double>& coord11, const std::vector<double>& coord12,
            const std::vector<double>& coord21, const std::vector<double>& coord22);
    
    ///Getters for constants
    double GetStretchingConstant(){return _kStretch;}
    double GetEqLength(){return _eqLength;}
    
    ///Setter and getter for parent linker
    void setLinker(Linker* linker) {_pLinker = linker;}
    Linker* getLinker() {return _pLinker;}
    
private:
    double _eqLength;
    double _kStretch;
    
    Linker* _pLinker; ///< ptr to parent linker
    
};

#endif /* defined(__CytoMech__MLinker__) */
