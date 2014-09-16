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
#include "FilamentDB.h"
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
    
    ///Print basic information about filaments
    void printFilaments() {
        
        for(auto &filament : *FilamentDB::Instance(FilamentDBKey())) {
            
            for (auto cylinder : filament->getCylinderVector()){
                
                auto x = cylinder->getMCylinder()->GetFirstBead()->coordinate;
                _outputFile<<x[0]<<","<<x[1]<<","<<x[2]<<"  ";
                
            }
            auto x = filament->getLastCylinder()->getMCylinder()->GetFirstBead()->coordinate;
            _outputFile<<x[0]<<","<<x[1]<<","<<x[2]<<std::endl;
        }
        _outputFile <<std::endl;
    }

};


#endif /* defined(__Cyto__Output__) */
