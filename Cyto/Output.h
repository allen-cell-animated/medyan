//
//  Output.h
//  Cyto
//
//  Created by James Komianos on 9/16/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Output__
#define __Cyto__Output__

#include <iostream>
#include <fstream>
#include "common.h"
#include "FilamentDB.h"
#include "CylinderDB.h"
#include "Bead.h"

class Output {

private:
    std::ofstream _outputFile; ///< output file being used
    
public:
    Output(std::string outputFileName) {
        _outputFile.open(outputFileName);
        if(!_outputFile.is_open()) {
            std::cout << "There was an error opening file " << outputFileName << " for output. Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Opening file " << outputFileName << std::endl;
    }
    ~Output() {_outputFile.close();}
    
    ///Print basic information about filaments (OLD)
    void printBasicSnapshot(int step);
    
    ///NEW OUTPUT
    void printSnapshot(int step);
    ///Print birth times of beads for each filament
    void printBirthTimes(int step);
    ///Print forces on beads for each filament
    void printForces(int step);
    ///Print stresses on beads for each filament
    void printStresses(int step);
};


#endif /* defined(__Cyto__Output__) */
